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
  t_float ratio;
  double phase;
} t_boid;

typedef struct _boids1 {
  t_object x_obj;

  int x_num_boids;
  t_boid *boids;

  // double x_phase;
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
  t_sample *in1 = (t_sample *)(w[2]); // fix the type
  t_sample *out = (t_sample *)(w[3]); // fix the type
  int n = (int)(w[4]);

  float *tab = cos_table;
  t_float conv = x->x_conv;
  double phase_0 = x->boids[0].phase;
  double phase_1 = x->boids[1].phase;
  double phase_2 = x->boids[2].phase;
  double phase_3 = x->boids[3].phase;

  if (!tab) return (w+5);

  while (n--) {
    t_sample root_f = *in1++;
    t_sample output = (t_sample)0.0;

    phase_0 += root_f * conv;
    double curphase_0 = phase_0;
    unsigned int idx_0 = (unsigned int)curphase_0;
    t_sample frac_0 = (t_sample)(curphase_0 - idx_0);
    idx_0 &= (WAVETABLE_SIZE - 1);
    t_sample f1_0 = tab[idx_0];
    t_sample f2_0 = tab[(idx_0 + 1) & (WAVETABLE_SIZE - 1)];
    output += f1_0 + frac_0 * (f2_0 - f1_0);

    phase_1 += root_f * x->boids[1].ratio * conv;
    double curphase_1 = phase_1;
    unsigned int idx_1 = (unsigned int)curphase_1;
    t_sample frac_1 = (t_sample)(curphase_1 - idx_1);
    idx_1 &= (WAVETABLE_SIZE - 1);
    t_sample f1_1 = tab[idx_1];
    t_sample f2_1 = tab[(idx_1 + 1) & (WAVETABLE_SIZE - 1)];
    output += f1_1 + frac_1 * (f2_1 - f1_1);

    phase_2 += root_f * x->boids[2].ratio * conv;
    double curphase_2 = phase_2;
    unsigned int idx_2 = (unsigned int)curphase_2;
    t_sample frac_2 = (t_sample)(curphase_2 - idx_2);
    idx_2 &= (WAVETABLE_SIZE - 1);
    t_sample f1_2 = tab[idx_2];
    t_sample f2_2 = tab[(idx_2 + 1) & (WAVETABLE_SIZE - 1)];
    output += f1_2 + frac_2 * (f2_2 - f1_2);

    phase_3 += root_f * x->boids[3].ratio * conv;
    double curphase_3 = phase_3;
    unsigned int idx_3 = (unsigned int)curphase_3;
    t_sample frac_3 = (t_sample)(curphase_3 - idx_3);
    idx_3 &= (WAVETABLE_SIZE - 1);
    t_sample f1_3 = tab[idx_3];
    t_sample f2_3 = tab[(idx_3 + 1) & (WAVETABLE_SIZE - 1)];
    output += f1_3 + frac_3 * (f2_3 - f1_3);

    *out++ = output * 0.25f;
  }

  while (phase_0 >= WAVETABLE_SIZE) phase_0 -= WAVETABLE_SIZE;
  while (phase_0 < 0) phase_0 += WAVETABLE_SIZE;
  x->boids[0].phase = phase_0;

  while (phase_1 >= WAVETABLE_SIZE) phase_1 -= WAVETABLE_SIZE;
  while (phase_1 < 0) phase_1 += WAVETABLE_SIZE;
  x->boids[1].phase = phase_1;

  while (phase_2 >= WAVETABLE_SIZE) phase_2 -= WAVETABLE_SIZE;
  while (phase_2 < 0) phase_2 += WAVETABLE_SIZE;
  x->boids[2].phase = phase_2;

  while (phase_3 >= WAVETABLE_SIZE) phase_3 -= WAVETABLE_SIZE;
  while (phase_3 < 0) phase_3 += WAVETABLE_SIZE;
  x->boids[3].phase = phase_3;

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
    // TODO: handle this better
    pd_error(x, "boids1~: failed to allocate memory to boids");
  }

  // srand(time(NULL)); // do this somewhere else (only call once)
  
  t_float r = 1.123f;
  for (int i = 0; i < x->x_num_boids; i++) {
    x->boids[i].ratio = r * i;
    x->boids[i].phase = (double)0.0;
    post("initialized boid[%d] with ratio %f", i, x->boids[i].ratio);
  }
}

static void *boids1_new(t_floatarg root_freq, t_floatarg num_boids)
{
  t_boids1 *x = (t_boids1 *)pd_new(boids1_class);

  // x->x_phase = (double)0.0;
  x->x_f = root_freq > 0 ? (t_float)root_freq : (t_float)220.0;
  // x->x_num_boids = num_boids > 0 ? (int)num_boids : 4;
  x->x_num_boids = 4;

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
                               A_DEFFLOAT, 0);

  class_addmethod(boids1_class, (t_method)boids1_dsp, gensym("dsp"), A_CANT, 0);
  CLASS_MAINSIGNALIN(boids1_class, t_boids1, x_f);
}


