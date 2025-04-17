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

typedef struct _boid {
  t_float ratio; // sets the ratio of a boid's frequency from the master
  // frequency.
  double wavetable_phase;
  t_float grain_duration_ms;
  t_float window_phase_inc;
  double window_phase;
  t_float cycle_start_pos;
  int active;
} t_boid;

typedef struct _boids3 {
  t_object x_obj;

  t_float x_sr;

  int x_num_boids;
  t_boid *boids;
  int x_num_active_boids;

  int x_cycle_ms;
  t_float x_cycle_pos;
  t_float x_ms_per_block;

  t_float x_conv;

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
    t_float window_samples = boid->grain_duration_ms * x->x_sr * (t_float)0.001;
    boid->window_phase_inc = (t_float)((t_float)WINDOWTABLE_SIZE / window_samples);
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

  t_float prev_cycle_pos = x->x_cycle_pos;
  t_float new_cycle_pos = prev_cycle_pos + x->x_ms_per_block;
  // rounding errors are fixed once per cycle; this may be good enough
  if (new_cycle_pos >= x->x_cycle_ms) new_cycle_pos = 0.0f;
  x->x_cycle_pos = new_cycle_pos;
  t_float next_cycle_pos = new_cycle_pos + x->x_ms_per_block;
  int active_count = 0;

  for (int i = 0; i < x->x_num_boids; i++) {
    t_boid *boid = &x->boids[i];
    if (!boid->active &&
        boid->cycle_start_pos > prev_cycle_pos &&
        boid->cycle_start_pos < next_cycle_pos) {
      boid->active = 1;
      x->x_num_active_boids++;
    }
    if (boid->active) active_count++;
  }

  t_float window_amp_scale = (t_float)((t_float)1.0 / (t_float)active_count); // tmp solution

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

        // todo: add boid amplitude attribute
        output += (f1 + frac * (f2 - f1)) * gsample * window_amp_scale;

        while (wphase >= WAVETABLE_SIZE) wphase -= WAVETABLE_SIZE;
        boid->wavetable_phase = wphase;

        while (gphase >= WINDOWTABLE_SIZE) gphase -= WINDOWTABLE_SIZE;
        boid->window_phase = gphase;

        if (boid->window_phase >= WINDOWTABLE_SIZE - boid->window_phase_inc) {
          boid->active = 0;
        }
      }
    }

    *out++ = output;
  }

  return (w+5);
}

static void boids3_dsp(t_boids3 *x, t_signal **sp)
{
  // conversion factor for sample rate
  x->x_conv = (t_float)WAVETABLE_SIZE / sp[0]->s_sr;
  x->x_sr = (t_float)sp[0]->s_sr;
  x->x_ms_per_block = ((t_float)1000.0 / x->x_sr) * sp[0]->s_length;
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

  int divisions = 32;
  t_float division_ms = (t_float)x->x_cycle_ms / (t_float)divisions;

  // efficiency isn't a big deal here, but maybe use a pointer:
  // t_boid *boid = &x->boids[i];
  for (int i = 0; i < x->x_num_boids; i++) {
    // initialize with random ratios between 0.5 and 2.0
    x->boids[i].ratio = (t_float)0.5 + (t_float)1.5 * ((t_float)rand() / RAND_MAX);
    x->boids[i].grain_duration_ms = (t_float)100.0; // hardcode for now
    x->boids[i].window_phase = (double)0.0;
    x->boids[i].wavetable_phase = (double)0.0;
    // initialize, then call update function from dsp method
    x->boids[i].window_phase_inc = (t_float)0.0;

    int division  = rand() % divisions;
    x->boids[i].cycle_start_pos = division * division_ms;
    x->boids[i].active = 0;
    post("boids3~: (debug) initialized boid[%d] with ratio %f", i, x->boids[i].ratio);
  }
}

static void *boids3_new(t_floatarg root_freq, t_floatarg num_boids)
{
  t_boids3 *x = (t_boids3 *)pd_new(boids3_class);

  x->x_f = root_freq > 0 ? root_freq : (t_float)220.0;
  x->x_num_boids = num_boids > 0 ? (int)num_boids : 4;
  x->x_sr = (t_float)0.0;

  x->x_cycle_ms = 8000;
  x->x_cycle_pos = (t_float)0.0;
  x->x_ms_per_block = (t_float)0.0;
  x->x_num_active_boids = 0;

  x->x_freq_inlet = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
  pd_float((t_pd *)x->x_freq_inlet, x->x_f);
  x->x_outlet = outlet_new(&x->x_obj, &s_signal);

  wavetable_init();
  windowtable_init();
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


