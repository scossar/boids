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
} t_cycle_division;

typedef struct _boid {
  t_float ratio;
  double wavetable_phase;
  t_float window_duration_ms;
  double window_phase_inc;
  double window_phase;
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

static void boids_free(t_boids4 *x)
{
  if (x->x_frequency_inlet) {
    inlet_free(x->x_frequency_inlet);
  }

  if (x->x_outlet) {
    outlet_free(x->x_outlet);
  }

  wavetable_free();
  windowtable_free();
}

static t_int *boids4_perform(t_int *w)
{
  t_boids4 *x = (t_boids4 *)(w[1]);
  t_sample *in1 = (t_sample *)(w[2]);
  t_sample *out1 = (t_sample *)(w[3]);
  int n = (int)(w[4]);

  t_float *wave_tab = cos_table;
  t_float *window_tab = window_table;
  t_float conv = x->x_conv;

  while (n--) {
    t_sample f = *in1++;
    t_sample output = (t_sample)0.0;

    *out1++ = output;
  }

  return (w+5);
}

static void boids4_dsp(t_boids4 *x, t_signal **sp)
{
  x->x_sr = sp[0]->s_sr;
  x->x_conv = (t_float)WAVETABLE_SIZE / x->x_sr;
  x->x_cycle_samples = (int)(x->x_cycle_ms * x->x_sr * 0.001f);

  dsp_add(boids4_perform, 4, x, sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_length);
}

static void *boids4_new(t_floatarg root_freq, t_floatarg num_boids)
{
  t_boids4 *x = (t_boids4 *)pd_new(boids4_class);

  x->x_f = root_freq > 0 ? root_freq : (t_float)220.0;
  x->x_num_boids = num_boids > 0 ? (int)num_boids : 4;
  x->x_sr = (t_float)0.0;

  x->x_cycle_ms = 4000;
  x->x_cycle_samples = 0;

  // tmp
  x->x_num_cycle_divisions = 4;

  x->boids = NULL;
  x->cycle_divisions = NULL;
  x->x_conv = (t_float)0.0;
  x->x_sr = (t_float)0.0;

  x->x_frequency_inlet = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
  pd_float((t_pd *)x->x_frequency_inlet, x->x_f);
  x->x_outlet = outlet_new(&x->x_obj, &s_signal);

  wavetable_init();
  windowtable_init();

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
}
