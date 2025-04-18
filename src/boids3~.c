#include "m_pd.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>

/*
* notes:
* - tracking_boid_phases.md
* - boid_window.md
* - boid_cycle_position.md
* - typecasting_in_c.md  
* - initializing_boids_and_cycle_divisions.md
*
*/

#define WAVETABLE_SIZE 4096 // 2^12 seems good enough
#define WINDOWTABLE_SIZE 4096 // use a separate variable in case I want to
// change it

static t_class *boids3_class = NULL;
static t_float *cos_table = NULL;
static int wavetable_reference_count = 0; // track wave tables
static t_float *window_table = NULL;
static int windowtable_reference_count = 0; // track window tables

typedef struct _cycle_division {
  t_float start_pos_ms; // the division's position in the cycle
  // pointer to start of array of boid indices that are currently associated with the division
  int *boid_indices;
  int num_boids; // number of boids associated with the division
  t_float amp_scale;
} t_cycle_division;

typedef struct _boid {
  // sets the ratio of the boid's frequency to the master (root) frequency
  t_float ratio;
  // NOTE: messing around:
  t_cycle_division *div;
  double wavetable_phase; // boid's current position in the wave table
  t_float window_duration_ms; // ms of the boid's (Hann) window
  // how much to increment the phase for each DSP step
  t_float window_phase_inc;
  double window_phase; // boid's current position in the window table
  // t_float cycle_start_pos;
  t_float window_deactivation_threshold;
  int active; // boolean (1 indicates active)
} t_boid;

typedef struct _boids3 {
  t_object x_obj;

  t_float x_sr; // system sample rate
  t_float x_ms_per_block; // ms per Pure Data block size

  int x_num_boids; // number of boids
  t_boid *boids; // pointer to start of boids array

  int x_cycle_ms; // length of cycle in ms
  int x_num_cycle_divisions; // number of divisions in the cycle
  t_cycle_division *cycle_divisions; // pointer to start of divisions array

  t_float x_cycle_pos_ms; // current possition in the cycle

  t_float x_conv; // frequency convesion for system sample rate

  t_inlet *x_freq_inlet;
  t_outlet *x_outlet;
  t_float x_f;
} t_boids3;

static void wavetable_init(void)
{
  if (cos_table == NULL) {
    // adding an extra sample for 2 point interpolation
    cos_table = (t_float *)getbytes(sizeof(t_float) * (WAVETABLE_SIZE + 1));
    if (cos_table) {
      // <= accounts for the extra sample here
      for (int i = 0; i <= WAVETABLE_SIZE; i++) {
        cos_table[i] = (t_float)(cos((i * (t_float)2.0 * (t_float)M_PI) / (t_float)WAVETABLE_SIZE));
      }
      post("boids3~: initialized cosine table of size %d (plus 1 sample)", WAVETABLE_SIZE);
    } else {
      post("boids3~ error: failed to allocate memory for cosine table");
      return;
    }
  }
  wavetable_reference_count++;
}

// start with just a Hann window
static void windowtable_init(void)
{
  if (window_table == NULL) {
    // extra sample for 2 point interpolation
    window_table = (t_float *)getbytes(sizeof(t_float) * (WINDOWTABLE_SIZE + 1));
    if (window_table) {
      for (int i = 0; i <= WINDOWTABLE_SIZE; i++) {
        window_table[i] = (t_float)(0.5 *
          (1.0 - cos((i * 2.0 * (t_float)M_PI) / (t_float)WINDOWTABLE_SIZE)));
      }
      post("boids3~: initialized window table of size %d (plus 1 sample)", WINDOWTABLE_SIZE);
    } else {
      post("boids3~: error: failed to allocate memory for window table");
      return;
    }
  }
  windowtable_reference_count++;
}

static void update_boid_parameters(t_boids3 *x)
{
  for (int i = 0; i < x->x_num_boids; i++) {
    t_boid *boid = &x->boids[i]; // use a pointer to avoid copying
    t_float window_samples = boid->window_duration_ms * x->x_sr * (t_float)0.001;
    boid->window_phase_inc = (t_float)((t_float)WINDOWTABLE_SIZE / window_samples);
    boid->window_deactivation_threshold = (t_float)WINDOWTABLE_SIZE - boid->window_phase_inc;
  }
}

static void wavetable_free(void)
{
  wavetable_reference_count--;
  if (wavetable_reference_count <= 0 && cos_table != NULL) {
    // NOTE: extra sample used for interpolation
    freebytes(cos_table, sizeof(t_float) * (WAVETABLE_SIZE + 1));
    cos_table = NULL;
    post("boids3~: freed cosine table");
    wavetable_reference_count = 0; // just to be safe
  }
}

static void windowtable_free(void)
{
  windowtable_reference_count--;
  if (windowtable_reference_count <= 0 && window_table != NULL) {
    // NOTE: extra sample used for interpolation
    freebytes(window_table, sizeof(t_float) * (WINDOWTABLE_SIZE + 1));
    window_table = NULL;
    post("boids3~: freed window table");
    windowtable_reference_count = 0;
  }
}

// add a boid to a division
static void add_boid_to_division(t_cycle_division *div, int boid_idx) {
  div->boid_indices[div->num_boids] = boid_idx;
  div->num_boids++;
  div->amp_scale++;
}

static void remove_boid_from_division(t_cycle_division *div, int boid_idx) {
  // just in case
  if (div->num_boids <= 0) return;

  // find the boid
  for (int i = 0; i < div->num_boids; i++) {
    if (div->boid_indices[i] == boid_idx) {
      // move last boid to this position, unless it's already the last boid
      if (i < div->num_boids - 1) {
        div->boid_indices[i] = div->boid_indices[div->num_boids - 1];
      }
      div->num_boids--;
      if (div->amp_scale > 1) div->amp_scale--;
      return;
    }
  }
}

static void activate_division_boids(t_boids3 *x, t_cycle_division *div)
{
  t_boid *boids = x->boids;
  int *indices = div->boid_indices;
  int count = div->num_boids;

  for (int i = 0; i < count; i++) {
    boids[indices[i]].active = 1;
  }
}

static void deactivate_division_boids(t_boids3 *x, t_cycle_division *div)
{
  t_boid *boids = x->boids;
  int *indices = div->boid_indices;
  int count = div->num_boids;

  for (int i = 0; i < count; i++) {
    boids[indices[i]].active = 0;
  }
}

static t_int *boids3_perform(t_int *w)
{
  t_boids3 *x = (t_boids3 *)(w[1]);
  t_sample *in1 = (t_sample *)(w[2]);
  t_sample *out = (t_sample *)(w[3]);
  int n = (int)(w[4]);

  t_float *tab = cos_table;
  t_float *gtable = window_table;
  t_float conv = x->x_conv;
  int wmask = WAVETABLE_SIZE - 1;
  int gmask = WINDOWTABLE_SIZE - 1;

  t_float prev_cycle_pos = x->x_cycle_pos_ms;
  t_float new_cycle_pos = prev_cycle_pos + x->x_ms_per_block;


  if (new_cycle_pos >= x->x_cycle_ms) {
    x->x_cycle_pos_ms = (t_float)0.0;
    t_cycle_division *div = &x->cycle_divisions[0];
    activate_division_boids(x, div);
  } else {
    x->x_cycle_pos_ms = new_cycle_pos;
    t_float next_cycle_pos = new_cycle_pos + x->x_ms_per_block;
    for (int i = 0; i < x->x_num_cycle_divisions; i++) {
      t_cycle_division *div = &x->cycle_divisions[i];
      if (div->start_pos_ms > prev_cycle_pos && div->start_pos_ms < next_cycle_pos) {
        activate_division_boids(x, div);
      }
    }
  }

  if (!tab) return (w+5);

  while (n--) {
    t_sample root_f = *in1++;
    t_sample output = (t_sample)0.0;

    for (int i = 0; i < x->x_num_boids; i++) {
      t_boid *boid = &x->boids[i]; // use a pointer to avoid copying

      if (boid->active) {
        double wphase = boid->wavetable_phase;
        unsigned int idx = (unsigned int)wphase & wmask;
        t_sample frac = (t_sample)(wphase - idx);
        t_sample f1 = tab[idx];
        t_sample f2 = tab[idx + 1];
        wphase += boid->ratio * root_f * conv;

        double gphase = boid->window_phase;
        unsigned int gidx = (unsigned int)gphase & gmask;
        t_sample gfrac = (t_sample)(gphase - gidx);
        t_sample g1 = gtable[gidx];
        t_sample g2 = gtable[gidx + 1];
        t_sample gsample = g1 + gfrac * (g2 - g1);
        gphase += boid->window_phase_inc;
        t_sample scale = (t_sample)1.0 / (t_sample)(boid->div->amp_scale);

        output += (f1 + frac * (f2 - f1)) * gsample * scale;

        while (wphase >= WAVETABLE_SIZE) wphase -= WAVETABLE_SIZE;
        boid->wavetable_phase = wphase;

        while (gphase >= WINDOWTABLE_SIZE) gphase -= WINDOWTABLE_SIZE;
        boid->window_phase = gphase;

        // deactivating boids is happening at sample rate
        // activating boids is happening at control rate
        // this allows for boids to extend past the duration of a single cycle
        // division (not currently being done though)
        if (boid->window_phase >= boid->window_deactivation_threshold) {
          boid->active = 0;
          boid->window_phase = (t_float)0.0;
        }
      }
    }

    *out++ = output;
  }

  return (w+5);
}

// called once when dsp is enabled on Pure Data
static void boids3_dsp(t_boids3 *x, t_signal **sp)
{
  // conversion factor for sample rate
  x->x_conv = (t_float)WAVETABLE_SIZE / sp[0]->s_sr;
  x->x_sr = (t_float)sp[0]->s_sr;
  // ms per sample * Pd sample block size
  x->x_ms_per_block = ((t_float)1000.0 / x->x_sr) * sp[0]->s_length;
  // use system sample rate to set boid window_phase_inc and
  // window_deactivation_threshold to correct values
  update_boid_parameters(x);
  dsp_add(boids3_perform, 4, x, sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_length);
}

static void boids_init(t_boids3 *x)
{
  x->boids = (t_boid *)getbytes(sizeof(t_boid) * x->x_num_boids);
  if (x->boids == NULL) {
    pd_error(x, "boids3~: failed to allocate memory to boids");
    return;
  }

  for (int i = 0; i < x->x_num_boids; i++) {
    t_boid *boid = &x->boids[i];
    boid->ratio = (t_float)0.5 + (t_float)1.5 * ((t_float)rand() / RAND_MAX);
    boid->window_duration_ms = (t_float)50.0; // hardcoded for now
    boid->window_phase = (double)0.0;
    boid->wavetable_phase = (double)0.0;
    boid->window_phase_inc = (t_float)0.0;
    boid->window_deactivation_threshold = (t_float)0.0;
    boid->div = NULL; // maybe safer?
  }

  for (int i = 0; i < x->x_num_boids; i++) {
    t_boid *boid = &x->boids[i];
    int division_idx = rand() % x->x_num_cycle_divisions;
    t_cycle_division *div = &x->cycle_divisions[division_idx];
    // add the boid's index to the indices array at position num_boids (next
    // available slot)
    div->boid_indices[div->num_boids] = i;
    boid->div = div; // experiment (associate by pointer vs index)
    div->num_boids++;
    div->amp_scale++;
  }
}

// TODO: return int or boolean values to indicate success/failure
static void cycle_divisions_init(t_boids3 *x)
{
  x->cycle_divisions = (t_cycle_division *)getbytes(sizeof(t_cycle_division)
                                                          * x->x_num_cycle_divisions);
  if (x->cycle_divisions == NULL) {
    pd_error(x, "boids3~: failed to allocate memory to cycle_divisions");
    return;
  }

  t_float division_ms = (t_float)x->x_cycle_ms / (t_float)x->x_num_cycle_divisions;
  for (int i = 0; i < x->x_num_cycle_divisions; i++){
    t_cycle_division *div = &x->cycle_divisions[i];
    // getbytes calls calloc, so boid_indices are being initialized to 0
    // since 0 is a valid boid indice, it would probably be better to initialize
    // to -1
    div->boid_indices = (int *)getbytes(sizeof(int) * x->x_num_boids);
    if (div->boid_indices == NULL) {
      pd_error(x, "boids3~: failed to allocate memory for cycle_division boid_indices");
      return;
    }
    div->num_boids = 0;
    div->amp_scale = (t_float)1.0;
    div->start_pos_ms = division_ms * (t_float)i;
  }
}

static void *boids3_new(t_floatarg root_freq, t_floatarg num_boids)
{
  t_boids3 *x = (t_boids3 *)pd_new(boids3_class);

  x->x_f = root_freq > 0 ? root_freq : (t_float)220.0;
  x->x_num_boids = num_boids > 0 ? (int)num_boids : 4;
  x->x_sr = (t_float)0.0;

  x->x_cycle_ms = 4000;
  x->x_ms_per_block = (t_float)0.0;

  // tmp
  x->x_num_cycle_divisions = 4;

  // initialize all pointers to NULL so that partial initialization failures can
  // be handled by the free function
  x->boids = NULL;
  x->cycle_divisions = NULL;

  x->x_freq_inlet = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
  pd_float((t_pd *)x->x_freq_inlet, x->x_f);
  x->x_outlet = outlet_new(&x->x_obj, &s_signal);

  // TODO: return ints on success/failure
  wavetable_init();
  windowtable_init();
  cycle_divisions_init(x);
  boids_init(x);

   static int seed_initialized = 0;
   if (!seed_initialized) {
     srand((unsigned int)time(NULL));
     seed_initialized = 1;
   }

  return (void *)x;
}

static void boids3_free(t_boids3 *x)
{
  if (x->x_freq_inlet) {
    inlet_free(x->x_freq_inlet);
  }

  if (x->x_outlet) {
    outlet_free(x->x_outlet);
  }

  if (x->cycle_divisions != NULL) {
    for (int i = 0; i < x->x_num_cycle_divisions; i++) {
      t_cycle_division *div = &x->cycle_divisions[i];
      if (div->boid_indices != NULL) {
        freebytes(div->boid_indices, sizeof(int) * x->x_num_boids);
        div->boid_indices = NULL;
      }
    }
    freebytes(x->cycle_divisions, sizeof(t_cycle_division) * x->x_num_cycle_divisions);
    x->cycle_divisions = NULL;
  }

  if (x->boids != NULL) {
    for (int i = 0; i < x->x_num_boids; i++) {
      if (x->boids[i].div != NULL) {
        freebytes(x->boids[i].div, sizeof(t_cycle_division));
        x->boids[i].div = NULL;
      }
    }
    freebytes(x->boids, sizeof(t_boid) * x->x_num_boids);
    x->boids = NULL;
  }

  // decrease reference counts and possibly free tables
  wavetable_free();
  windowtable_free();
}

void boids3_tilde_setup(void)
{
  boids3_class = class_new(gensym("boids3~"),
                               (t_newmethod)boids3_new,
                               (t_method)boids3_free,
                               sizeof(t_boids3),
                               CLASS_DEFAULT,
                               A_DEFFLOAT, A_DEFFLOAT, 0);

  class_addmethod(boids3_class, (t_method)boids3_dsp, gensym("dsp"), A_CANT, 0);
  CLASS_MAINSIGNALIN(boids3_class, t_boids3, x_f);
}


