#include "m_pd.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>

// #define WAVETABLE_SIZE 16384 // 2^14
#define WAVETABLE_SIZE 4096 // 2^12 might be good enough

static t_class *boids1_class = NULL;
static float *cos_table = NULL; // shared wavetable
static int table_reference_count = 0; // track how many instances exist

typedef struct _boid {
  t_float ratio; // sets the ratio of a boid's frequency from the master
  // frequency.
  double phase;
} t_boid;

typedef struct _boids1 {
  t_object x_obj;

  int x_num_boids;
  t_boid *boids;

  t_float x_conv;
  t_inlet *x_freq_inlet;
  t_outlet *x_outlet;
  t_float x_f;
} t_boids1;

static void wavetable_init(void)
{
  if (cos_table == NULL) {
    cos_table = (float *)getbytes(sizeof(float) * (WAVETABLE_SIZE ));
    if (cos_table) {
      for (int i = 0; i < WAVETABLE_SIZE; i++) {
        cos_table[i] = cosf((i * 2.0f * (float)M_PI) / (float)WAVETABLE_SIZE);
      }
      post("boids1~: initialized cosine table of size %d", WAVETABLE_SIZE);
    } else {
      post("boids1~ error: failed to allocate memory for cosine table");
    }
  }
  table_reference_count++;
}

static void wavetable_free(void)
{
  table_reference_count--;
  if (table_reference_count <= 0 && cos_table != NULL) {
    freebytes(cos_table, sizeof(float) * (WAVETABLE_SIZE));
    cos_table = NULL;
    post("boids1~: freed cosine table");
    table_reference_count = 0; // just to be safe
  }
}

static t_int *boids1_perform(t_int *w)
{
  t_boids1 *x = (t_boids1 *)(w[1]);
  t_sample *in1 = (t_sample *)(w[2]);
  t_sample *out = (t_sample *)(w[3]);
  int n = (int)(w[4]);

  float *tab = cos_table;
  t_float conv = x->x_conv;
  int mask = WAVETABLE_SIZE - 1;
  t_float amp_scale = 1.0f / x->x_num_boids; // tmp solution

  if (!tab) return (w+5);

  while (n--) {
    t_sample root_f = *in1++;
    t_sample output = (t_sample)0.0;

    for (int i = 0; i < x->x_num_boids; i++) {
      t_boid boid = x->boids[i];
      double curphase = boid.phase;
      unsigned int idx = (unsigned int)curphase;
      t_sample frac = (t_sample)(curphase - idx);
      idx &= mask; // keep the index within the wavetable
      t_sample f1 = tab[idx];
      t_sample f2 = tab[(idx + 1) & mask]; // in case idx + 1 is out of the
      // wavetable
      output += f1 + frac * (f2 - f1); // add a boid amplitude attribute
      curphase += boid.ratio * root_f * conv;
      x->boids[i].phase = curphase; // can't use the bit mask on float values
    }

    *out++ = output * amp_scale;
  }

  for (int i = 0; i < x->x_num_boids; i++) {
    double phase = x->boids[i].phase;
    while (phase >= WAVETABLE_SIZE) phase -= WAVETABLE_SIZE;
    while (phase < 0) phase += WAVETABLE_SIZE;
    x->boids[i].phase = phase;
  }

  return (w + 5);
}

static void boids1_dsp(t_boids1 *x, t_signal **sp)
{
  // calculate the conversion factor for this sample rate
  x->x_conv = (float)WAVETABLE_SIZE / sp[0]->s_sr;

  dsp_add(boids1_perform, 4, x, sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_length);
}

static void boids_init(t_boids1 *x)
{
  x->boids = (t_boid *)getbytes(sizeof(t_boid) * x->x_num_boids);
  if (x->boids == NULL) {
    pd_error(x, "boids1~: failed to allocate memory to boids");
  }

  // srand(time(NULL)); // do this somewhere else (only call once)
  
  t_float r = 1.123f;
  for (int i = 0; i < x->x_num_boids; i++) {
    x->boids[i].ratio = r * (i + 1);
    x->boids[i].phase = (double)0.0;
    post("initialized boid[%d] with ratio %f", i, x->boids[i].ratio);
  }
}

static void *boids1_new(t_floatarg root_freq, t_floatarg num_boids)
{
  t_boids1 *x = (t_boids1 *)pd_new(boids1_class);

  x->x_f = root_freq > 0 ? (t_float)root_freq : (t_float)220.0;
  x->x_num_boids = num_boids > 0 ? (int)num_boids : 4;

  x->x_freq_inlet = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
  pd_float((t_pd *)x->x_freq_inlet, x->x_f);
  x->x_outlet = outlet_new(&x->x_obj, &s_signal);

  wavetable_init();
  boids_init(x);

  return (void *)x;
}

static void boids1_free(t_boids1 *x)
{
  if (x->x_freq_inlet) {
    inlet_free(x->x_freq_inlet);
  }

  if (x->x_outlet) {
    outlet_free(x->x_outlet);
  }

  // decrease reference count and possibly free wavetable
  wavetable_free();
}

void boids1_tilde_setup(void)
{
  boids1_class = class_new(gensym("boids1~"),
                               (t_newmethod)boids1_new,
                               (t_method)boids1_free,
                               sizeof(t_boids1),
                               CLASS_DEFAULT,
                               A_DEFFLOAT, A_DEFFLOAT, 0);

  class_addmethod(boids1_class, (t_method)boids1_dsp, gensym("dsp"), A_CANT, 0);
  CLASS_MAINSIGNALIN(boids1_class, t_boids1, x_f);
}


