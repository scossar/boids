#include "m_pd.h"

typedef struct _cycle_division {
  t_float start_position; 
  int num_boids;
  int *boid_indices;
  int capacity; // max boids (used to allocate/dealocate memory)
} t_cycle_division;

typedef struct _cycle {
  int num_divisions;
  t_cycle_division *divisions;
} t_cycle;

typedef struct _boid {
  t_float ratio; // ratio of boid frequency to master (tonic) frequency
  double wavetable_phase; // current phase in wave table
  t_float grain_duration_ms; // duration of grain and (Hann) window
  t_float window_phase_inc; // how far window phase advances for each DSP step
  double window_phase; // current position in (Hann) window
  int active; // non-zero indicates active
} t_boid;

// NOTE: the 2 in boids2 is just a temporary hack to allow me to create boids1, boids2, boids3,
// ... otherwise it's meaningless
typedef struct _boids2 {
  t_object x_obj;
  
  t_sample x_sr; // system sample rate
  t_float x_ms_per_block; // ms per system block size for system sample rate

  int x_num_boids; // total number of boids
  t_boid *boids; // pointer to start of boids memory location

  // assuming (for now) there's a single cycle associated with the object
  t_cycle *cycle; // pointer to the start of the cycle memory location
  int x_cycle_ms; // number of milliseconds per cycle (ie, BPM)
  t_float x_cycle_position; // current position in cycle

  t_float x_conv; // converts a frequency into a phase increment

  t_inlet *x_freq_inlet;
  t_outlet *x_outlet;
  t_float x_f;
} t_boids2;
