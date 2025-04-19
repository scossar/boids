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
* - refactoring_boids_external.md
*
*/

#define WAVETABLE_SIZE 4096 // 2^12
#define WINDOWTABLE_SIZE 4096

static t_class *boids4_class = NULL;
static t_float *cos_table = NULL;
static int wavetable_reference_count = 0; // track wave tables
static t_float *window_table = NULL;
static int windowtable_reference_count = 0; // track window tables

typedef struct _cycle_division {
  int start_sample;
  int end_sample;
  int num_boids;
  int *boid_indices;
  struct _cycle_division *next_division;
} t_cycle_division;

typedef struct _boid {
  t_float ratio;
  double wavetable_phase;
  t_float window_duration_ms;
  double window_phase_inc;
  double window_phase;
  t_float window_deactivation_threshold;
  
  t_cycle_division *current_division;
  int active;
} t_boid;

typedef struct _boids4 {
  t_object x_obj;

  t_float x_sr; // system sample rate

  int x_num_boids;
  t_boid *boids;

  int x_cycle_ms; // for user input; will be converted to samples internally
  int x_cycle_samples; // total samples in cycle
  int x_cycle_phase; // current sample in cycle

  int x_num_cycle_divisions;
  t_cycle_division *cycle_divisions;

  int x_current_division_index;

  t_float x_f; // root frequency
  t_float x_conv; // frequency conversion for sample rate

  t_inlet *x_frequency_inlet;
  t_outlet *x_outlet;
} t_boids4;

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
      post("boids4~: initialized cosine table of size %d (plus 1 sample)", WAVETABLE_SIZE);
    } else {
      post("boids4~ error: failed to allocate memory for cosine table");
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
      post("boids4~: initialized window table of size %d (plus 1 sample)", WINDOWTABLE_SIZE);
    } else {
      post("boids4~: error: failed to allocate memory for window table");
      return;
    }
  }
  windowtable_reference_count++;
}

static void wavetable_free(void)
{
  wavetable_reference_count--;
  if (wavetable_reference_count <= 0 && cos_table != NULL) {
    // NOTE: extra sample used for interpolation
    freebytes(cos_table, sizeof(t_float) * (WAVETABLE_SIZE + 1));
    cos_table = NULL;
    post("boids4~: freed cosine table");
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
    post("boids4~: freed window table");
    windowtable_reference_count = 0;
  }
}

static void cycle_divisions_init(t_boids4 *x)
{
  x->cycle_divisions = (t_cycle_division *)getbytes(sizeof(t_cycle_division) *
                                                    x->x_num_cycle_divisions);
  if (x->cycle_divisions == NULL) {
    pd_error(x, "boids4~: failed to allocate memory for cycle divisions");
    return;
  }

  for (int i = 0; i < x->x_num_cycle_divisions; i++) {
    x->cycle_divisions[i].boid_indices = (int *)getbytes(sizeof(int) * x->x_num_boids);
    if (x->cycle_divisions[i].boid_indices == NULL) {
      pd_error(x, "boids4~: failed to allocate memory for cycle_division[%d].boid_indeces", i);
      return;
    }

    int next_idx = (i + 1) % x->x_num_cycle_divisions;
    x->cycle_divisions[i].next_division = &x->cycle_divisions[next_idx];
  }
}

// called from DSP method
static void cycle_division_update(t_boids4 *x) {
  // the values are ints, but just in case:
  int samples_per_division = (int)x->x_cycle_samples / (int)x->x_num_cycle_divisions;
  int remainder = (int)x->x_cycle_samples % (int)x->x_num_cycle_divisions;

  int current_start = 0;
  for (int i = 0; i < x->x_num_cycle_divisions; i++) {
    x->cycle_divisions[i].start_sample = current_start;

    // deal with possible remainder (largest possible remainder is num_divisions
    // - 1)
    int this_division_samples = samples_per_division;
    if (i < remainder) {
      this_division_samples++;
    }
    current_start += this_division_samples;
    // (last division's end sample will be x_cycle_samples - 1)
    x->cycle_divisions[i].end_sample = current_start - 1;
  }
}

static void boids_init(t_boids4 *x)
{
  x->boids = (t_boid *)getbytes(sizeof(t_boid) * x->x_num_boids);
  if (x->boids == NULL) {
    pd_error(x, "boids4~: failed to allocate memory for boids");
    return;
  }

  for (int i = 0; i < x->x_num_boids; i++) {
    t_boid *boid = &x->boids[i];
    boid->ratio = (t_float)0.5 + (t_float)4.5 * ((t_float)rand() / RAND_MAX);
    boid->window_duration_ms = (t_float)220.0;
    boid->window_phase = (double)0.0;
    boid->window_phase_inc = (double)0.0;
    boid->wavetable_phase = (double)0.0;
    boid->active = 0;

    int division_idx = rand() % x->x_num_cycle_divisions;
    t_cycle_division *div = &x->cycle_divisions[division_idx];

    // add boid idx to division
    div->boid_indices[div->num_boids] = i;
    div->num_boids++;

    // set bidirectional reference
    boid->current_division = div;
  }
}

static void boids_update(t_boids4 *x)
{
  for (int i = 0; i < x->x_num_boids; i++) {
    t_boid *boid = &x->boids[i];
    t_float window_samples = boid->window_duration_ms * x->x_sr * (t_float)0.001;
    boid->window_phase_inc = (t_float)((t_float)WINDOWTABLE_SIZE / window_samples);
    boid->window_deactivation_threshold = (t_float)WINDOWTABLE_SIZE - boid->window_phase_inc;
  }
}

static void activate_division_boids(t_boids4 *x, t_cycle_division *div)
{
  t_boid *boids = x->boids;
  int *indices = div->boid_indices;
  int count = div->num_boids;

  for (int i = 0; i < count; i++) {
    boids[indices[i]].active = 1;
  }
}

static void move_boid_to_division(t_boids4 *x, int boid_idx, t_cycle_division *new_div)
{
  t_boid *boid = &x->boids[boid_idx];
  t_cycle_division *old_div = boid->current_division;
  // t_cycle_division *new_div = &x->cycle_divisions[new_div_idx];

  int found = 0;
  for (int i = 0; i < old_div->num_boids; i++) {
    if (found) {
      old_div->boid_indices[i-1] = old_div->boid_indices[i]; 
    } else if (old_div->boid_indices[i] == boid_idx) {
      found = 1;
    }
  }

  if (found) {
    old_div->num_boids--;
    new_div->boid_indices[new_div->num_boids] = boid_idx;
    new_div->num_boids++;

    boid->current_division = new_div;
    boid->active = 0; // on the assumption this is called before playing current
    // div
  }
}

// note: the move function can be changed to accept two division pointers (old
// and new)
static void apply_boid_threshold_rule(t_boids4 *x, t_cycle_division *div, int threshold)
{
  if (div->num_boids > threshold) {
    int random_idx = rand() % div->num_boids;
    int boid_idx = div->boid_indices[random_idx];
    move_boid_to_division(x, boid_idx, div->next_division);
  }
}

static void apply_div_harmonizing_rule(t_boids4 *x, t_cycle_division *div)
{
  for (int i = 0; i < div->num_boids; i++) {
    t_boid *boid = &x->boids[div->boid_indices[i]]; 
    t_float ratio = boid->ratio;
    t_float frac = ratio - (int)ratio;
    t_float scale = frac * (t_float)0.01;
    scale = (frac > (t_float)0.5) ? scale : -scale;
    boid->ratio += scale;
  }
}

static t_int *boids4_perform(t_int *w)
{
  t_boids4 *x = (t_boids4 *)(w[1]);
  t_sample *in1 = (t_sample *)(w[2]);
  t_sample *out1 = (t_sample *)(w[3]);
  int n = (int)(w[4]);

  t_float *wave_tab = cos_table;
  t_float *window_tab = window_table;
  int wave_mask = WAVETABLE_SIZE - 1;
  int window_mask = WINDOWTABLE_SIZE - 1;
  t_float conv = x->x_conv;

  int current_phase = x->x_cycle_phase;
  int next_phase = current_phase + n;
  // if (next_phase >= x->x_cycle_samples) next_phase -= x->x_cycle_samples;
  if (next_phase >= x->x_cycle_samples) next_phase = 0;

  int current_division_index = x->x_current_division_index;
  t_cycle_division *current_division = &x->cycle_divisions[current_division_index];
  int current_division_end = current_division->end_sample;

  if (current_division_end - current_phase < n) {
    int next_division_index = current_division_index + 1;
    if (next_division_index >= x->x_num_cycle_divisions) next_division_index = 0;
    x->x_current_division_index = next_division_index;
    activate_division_boids(x, &x->cycle_divisions[next_division_index]);
    apply_boid_threshold_rule(x, &x->cycle_divisions[next_division_index], 1);
    apply_div_harmonizing_rule(x, &x->cycle_divisions[next_division_index]);
  } else if (current_phase == 0 && current_division_index == 0) {
    activate_division_boids(x, &x->cycle_divisions[0]);
    apply_boid_threshold_rule(x, &x->cycle_divisions[0], 1);
    apply_div_harmonizing_rule(x, &x->cycle_divisions[0]);
  }

  // get a pointer to the current cycle_division
  t_cycle_division *div = &x->cycle_divisions[x->x_current_division_index];
  int num_boids = div->num_boids;
  t_float amp_scale = (num_boids > 1) ? (t_float)1.0 / (t_float)num_boids : (t_float)1.0;

  while (n--) {
    t_sample f = *in1++;
    t_sample output = (t_sample)0.0;

    for (int i = 0; i < x->x_num_boids; i++) {
      t_boid *boid = &x->boids[i];

      if (boid->active) {
        double w_phase = boid->wavetable_phase;
        unsigned int w_idx = (unsigned int)w_phase & wave_mask;
        t_sample w_frac = (t_sample)(w_phase - (t_sample)w_idx);
        t_sample f1 = wave_tab[w_idx];
        t_sample f2 = wave_tab[w_idx + 1];
        w_phase += boid->ratio * f * conv;

        double g_phase = boid->window_phase;
        unsigned int g_idx = (unsigned int)g_phase & window_mask;
        t_sample g_frac = (t_sample)(g_phase - (t_sample)g_idx);
        t_sample g1 = window_tab[g_idx];
        t_sample g2 = window_tab[g_idx + 1];
        t_sample g_sample = g1 + g_frac * (g2 - g1);
        g_phase += boid->window_phase_inc;

        output += (f1 + w_frac * (f2 - f1)) * g_sample * amp_scale;

        while (w_phase >= WAVETABLE_SIZE) w_phase -= WAVETABLE_SIZE;
        boid->wavetable_phase = w_phase;

        while (g_phase >= WINDOWTABLE_SIZE) g_phase -= WINDOWTABLE_SIZE;
        boid->window_phase = g_phase;

        if (boid->window_phase >= boid->window_deactivation_threshold) {
          boid->active = 0;
          boid->window_phase = (t_float)0.0;
        }
      }
    }

    *out1++ = output;
  }

  x->x_cycle_phase = next_phase;
  return (w+5);
}

static void boids4_dsp(t_boids4 *x, t_signal **sp)
{
  x->x_sr = sp[0]->s_sr;
  x->x_conv = (t_float)WAVETABLE_SIZE / x->x_sr;
  x->x_cycle_samples = (int)(x->x_cycle_ms * x->x_sr * 0.001f);
  cycle_division_update(x);
  boids_update(x);
  dsp_add(boids4_perform, 4, x, sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_length);
}

static void boids_free(t_boids4 *x)
{
  if (x->x_frequency_inlet != NULL) {
    inlet_free(x->x_frequency_inlet);
  }

  if (x->x_outlet != NULL) {
    outlet_free(x->x_outlet);
  }

  if (x->cycle_divisions != NULL) {
    for (int i = 0; i < x->x_num_cycle_divisions; i++) {
      freebytes(x->cycle_divisions[i].boid_indices, sizeof(int) * x->x_num_boids);
      x->cycle_divisions[i].boid_indices = NULL;
    }

    freebytes(x->cycle_divisions, sizeof(t_cycle_division) * x->x_num_cycle_divisions);
    x->cycle_divisions = NULL;
  }

  if (x->boids != NULL) {
    for (int i = 0; i < x->x_num_boids; i++) {
      if (x->boids[i].current_division != NULL) {
        // the memory has already been freed (it's owned by cycle_divisions)
        // just set to NULL to be safe
        x->boids[i].current_division = NULL;
      }
    }

    freebytes(x->boids, sizeof(t_boid) * x->x_num_boids);
    x->boids = NULL;
  }

  wavetable_free();
  windowtable_free();
}

static void *boids4_new(t_floatarg root_freq, t_floatarg num_boids)
{
  t_boids4 *x = (t_boids4 *)pd_new(boids4_class);

  x->x_f = root_freq > 0 ? root_freq : (t_float)220.0;
  x->x_num_boids = num_boids > 0 ? (int)num_boids : 4;
  x->x_sr = (t_float)0.0;

  x->x_cycle_ms = 2000;
  x->x_cycle_samples = 0;

  // tmp
  x->x_num_cycle_divisions = 16;

  x->boids = NULL;
  x->cycle_divisions = NULL;
  x->x_conv = (t_float)0.0;
  x->x_sr = (t_float)0.0;
  x->x_current_division_index = 0;
  x->x_cycle_phase = 0;

  x->x_frequency_inlet = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
  pd_float((t_pd *)x->x_frequency_inlet, x->x_f);
  x->x_outlet = outlet_new(&x->x_obj, &s_signal);

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

void boids4_tilde_setup(void)
{
  boids4_class = class_new(gensym("boids4~"),
                           (t_newmethod)boids4_new,
                           (t_method)boids_free,
                           sizeof(t_boids4),
                           CLASS_DEFAULT,
                           A_DEFFLOAT, A_DEFFLOAT, 0);

  class_addmethod(boids4_class, (t_method)boids4_dsp, gensym("dsp"), A_CANT, 0);
  CLASS_MAINSIGNALIN(boids4_class, t_boids4, x_f);
}
