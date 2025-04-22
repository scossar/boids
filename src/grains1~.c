#include "m_pd.h"
#include <math.h>

/*
 * implements grain spreading feature
 * notes:
 * - interpolation_from_frac_calculation.md
 */

#define WAVETABLE_SIZE 4096 // 2^12
#define WINDOWTABLE_SIZE 4096
#define XTRASAMPS 1

static t_class *grains1_class = NULL;
static t_float *wavetable = NULL;
static int wavetable_reference_count = 0;
static t_float *windowtable = NULL;
static int windowtable_reference_count = 0;

typedef struct _grain {
  t_float ratio;
  double wavetable_phase;

  t_float window_duration_ms;
  t_float window_samples;
  double windowtable_phase_inc;
  double windowtable_phase;
  int active;
  int start_countdown_samps; // samples until grain should become active
  int start_countdown_samps_remaining;
} t_grain;

typedef struct _grains1 {
  t_object x_obj;

  t_float x_sr;
  int x_pd_block_size;

  t_float x_f; // root frequency
  t_float x_conv; // frequency conversion for sample rate

  t_float x_init_spread;

  t_grain *grains;
  int x_num_grains;
  t_float x_num_grains_factor; // (t_float)1.0 / (t_float)x->x_num_grains

  t_inlet *x_freq_inlet;
  t_inlet *x_spread_inlet;
  t_outlet *x_outlet;
} t_grains1;

static void wavetable_init(void)
{
  if (wavetable == NULL) {
    wavetable = (t_float *)getbytes(sizeof(t_float) * (WAVETABLE_SIZE + XTRASAMPS));
    if (wavetable) {
      for (int i = 0; i < WAVETABLE_SIZE + XTRASAMPS; i++) {
        wavetable[i] = (t_float)(cos((i * (t_float)2.0 * (t_float)M_PI) / (t_float)WAVETABLE_SIZE));
      }
      post("grains1~: initialized cosine table of size %d plus %d samples", WAVETABLE_SIZE, XTRASAMPS);

    } else {
      post("grains1~: failed to allocate memory for wavetable");
      return;
    }
  }
  wavetable_reference_count++;
}

// Hann window
static void windowtable_init(void)
{
  if (windowtable == NULL) {
    // extra sample for 2 point interpolation
    windowtable = (t_float *)getbytes(sizeof(t_float) * (WINDOWTABLE_SIZE + XTRASAMPS));
    if (windowtable) {
      for (int i = 0; i < WINDOWTABLE_SIZE + XTRASAMPS; i++) {
        windowtable[i] = (t_float)(0.5 *
          (1.0 - cos((i * 2.0 * (t_float)M_PI) / (t_float)WINDOWTABLE_SIZE)));
      }
      post("grains1~: initialized window table of size %d + %d samples", WINDOWTABLE_SIZE, XTRASAMPS);
    } else {
      post("grains1~: failed to allocate memory for windowtable");
      return;
    }
  }
  windowtable_reference_count++;
}

static void wavetable_free(void)
{
  wavetable_reference_count--;
  if (wavetable_reference_count <= 0 && wavetable != NULL) {
    freebytes(wavetable, sizeof(t_float) * (WAVETABLE_SIZE + XTRASAMPS));
    wavetable = NULL;
    post("grains1~: freed wavetable");
    wavetable_reference_count = 0; // just to be safe
  }
}

static void windowtable_free(void)
{
  windowtable_reference_count--;
  if (windowtable_reference_count <= 0 && windowtable != NULL) {
    freebytes(windowtable, sizeof(t_float) * (WINDOWTABLE_SIZE + XTRASAMPS));
    windowtable = NULL;
    post("grains1~: freed windowtable");
    windowtable_reference_count = 0;
  }
}

static void grains1_free(t_grains1 *x)
{
  if (x->grains) {
    freebytes(x->grains, sizeof(t_grain) * x->x_num_grains);
    x->grains = NULL;
  }

  wavetable_free();
  windowtable_free();
}

static void grains_init(t_grains1 *x)
{
  x->grains = (t_grain *)getbytes(sizeof(t_grain) * x->x_num_grains);
  if (x->grains == NULL) {
    pd_error(x, "grains1~: failed to allocate memory for grains");
    return;
  }

  // tmp hack
  int max_harmonics = 11;
  for (int i = 0; i < x->x_num_grains; i++) {
    t_grain *grain = &x->grains[i];
    grain->ratio = (t_float)((i+1) % max_harmonics);
    grain->window_duration_ms = (t_float)125.0; // tmp hardcoded
    grain->wavetable_phase = (double)0.0;
    grain->windowtable_phase = (double)0.0;
    grain->windowtable_phase_inc = (double)0.0;
    grain->start_countdown_samps = 0;
    grain->start_countdown_samps_remaining = 0;
    grain->active = 0;
  }
}

// called from DSP initializer
static void grains_update(t_grains1 *x, t_float spread_factor)
{
  spread_factor = (spread_factor <= 0) ? (t_float)1.0 : spread_factor;
  for (int i = 0; i < x->x_num_grains; i++) {
    t_grain *grain = &x->grains[i];
    t_float window_samples = grain->window_duration_ms * x->x_sr * (t_float)0.001;
    grain->window_samples = window_samples;
    grain->windowtable_phase_inc = (double)(t_float)WINDOWTABLE_SIZE / window_samples;
    // x_num_grains_factor is 1.0 / x_num_grains
    int spacing = grain->window_samples * x->x_num_grains_factor * spread_factor;
    grain->start_countdown_samps = i * spacing;
    grain->start_countdown_samps_remaining = grain->start_countdown_samps;
  }
}

// passing the t_grains1 object as an arg might be adding unnecessary overhead
// instead, x->x_num_grains_factor could be passed as an argument
static void grain_distribute(t_grains1 *x, t_grain *grain, int i, t_float spread_factor)
{
  spread_factor = (spread_factor <= 0) ? (t_float)1.0 : spread_factor;
  // x_num_grains_factor allows for multiplication instead of requiring division
  int spacing = grain->window_samples * x->x_num_grains_factor * spread_factor;
  // using (i + 1) so that x->grains[0] isn't excluded
  grain->start_countdown_samps = (i + 1) * spacing;
  // not sure that `delta` is an appropriate variable name here
  // the goal is to accunt for the updated start_countdown_samps value without
  // overwriting samples that have already been decremented
  // not sure that's necessary
  int delta = grain->start_countdown_samps - grain->start_countdown_samps_remaining;
  grain->start_countdown_samps_remaining += delta;
}

static t_int *grains1_perform(t_int *w)
{
  t_grains1 *x = (t_grains1 *)(w[1]);
  t_sample *in1 = (t_sample *)(w[2]);
  t_sample *in2 = (t_sample *)(w[3]);
  t_sample *out1 = (t_sample *)(w[4]);
  int n = (int)(w[5]);

  t_float *wave_t = wavetable;
  t_float *win_t = windowtable;
  int wave_mask = WAVETABLE_SIZE - 1;
  int win_mask = WINDOWTABLE_SIZE - 1;
  t_float conv = x->x_conv;

  // tmp solution
  t_sample grains_scale = (t_float)1.0 / (t_float)x->x_num_grains;

  while (n--) {
    t_sample f = *in1++;
    t_sample spread = *in2++;
    if (PD_BIGORSMALL(f)) f = 0;

    t_sample output = (t_sample)0.0;
    
    for (int i = 0; i < x->x_num_grains; i++) {
      t_grain *grain = &x->grains[i];

      if (!grain->active) {
        if (grain->start_countdown_samps_remaining <= 0) {
          grain->active = 1;
          grain_distribute(x, grain, i, spread);
        } else {
          grain->start_countdown_samps_remaining--;
          continue;
        }
      }

      double wave_phase = grain->wavetable_phase;
      while (wave_phase >= WAVETABLE_SIZE) wave_phase -= WAVETABLE_SIZE;
      unsigned int wave_idx = (unsigned int)wave_phase & wave_mask;
      t_sample wave_frac = (t_sample)(wave_phase - (t_sample)wave_idx);
      t_sample f1 = wave_t[wave_idx];
      t_sample f2 = wave_t[wave_idx+1];
      wave_phase += grain->ratio * f * conv;
      grain->wavetable_phase = wave_phase;

      double win_phase = grain->windowtable_phase;
      while (win_phase >= WINDOWTABLE_SIZE) win_phase -= WINDOWTABLE_SIZE;
      unsigned int win_idx = (unsigned int)win_phase & win_mask;
      t_sample win_frac = (t_sample)(win_phase - (t_sample)win_idx);
      t_sample win1 = win_t[win_idx];
      t_sample win2 = win_t[win_idx+1];
      t_sample window = win1 + win_frac * (win2 - win1);
      win_phase += grain->windowtable_phase_inc;
      grain->windowtable_phase = win_phase;
      if (win_phase >= WINDOWTABLE_SIZE) {
        grain->active = 0;
      }

      output += (f1 + wave_frac * (f2 - f1)) * window * grains_scale;
    }
    *out1++ = output;
  }

  return (w+6);
}

static void grains1_dsp(t_grains1 *x, t_signal **sp)
{
  x->x_sr = sp[0]->s_sr;
  x->x_pd_block_size = sp[0]->s_length;
  x->x_conv = (t_float)WAVETABLE_SIZE / x->x_sr;
  grains_update(x, (t_float)1.5);
  dsp_add(grains1_perform, 5, x, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[0]->s_length);
}

static void *grains1_new(t_floatarg root_freq, t_floatarg num_grains)
{
  t_grains1 *x = (t_grains1 *)pd_new(grains1_class);
  x->x_f = root_freq > 0 ? root_freq : (t_float)220.0;
  x->x_num_grains = num_grains > 0 ? (int)num_grains : 4;
  x->x_num_grains_factor = (t_float)1.0 / (t_float)x->x_num_grains;

  x->x_sr = (t_float)0.0;
  x->x_pd_block_size = 0;
  x->x_init_spread = (t_float)1.0;

  wavetable_init();
  windowtable_init();
  grains_init(x);

  x->x_freq_inlet = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
  pd_float((t_pd *)x->x_freq_inlet, x->x_f);
  x->x_spread_inlet = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
  pd_float((t_pd *)x->x_spread_inlet, x->x_init_spread);

  x->x_outlet = outlet_new(&x->x_obj, &s_signal);

  return (void *)x;
}

void grains1_tilde_setup(void)
{
  grains1_class = class_new(gensym("grains1~"),
                            (t_newmethod)grains1_new,
                            (t_method)grains1_free,
                            sizeof(t_grains1),
                            CLASS_DEFAULT,
                            A_DEFFLOAT, A_DEFFLOAT, 0);

  class_addmethod(grains1_class, (t_method)grains1_dsp, gensym("dsp"), A_CANT, 0);
  CLASS_MAINSIGNALIN(grains1_class, t_grains1, x_f);
}



