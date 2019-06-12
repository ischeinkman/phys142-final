#include "hermite_polynomial.h"
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <memory.h>
#include <stdlib.h>
#include <time.h>

#define min(a, b) ((a > b) ? b : a)

// Particle mass
#define m (1.0)

// HO Frequency
#define w (1.0)

// Reduced Planck's constant
#define HBAR (1.0)

// Boltzman Constant
#define k (1.0)

#define TEMP (HBAR/(k * TOTAL_TIME))
// number of discrete points of the polymer chain
#define N (500)

// number of runs of the monte carlo path simulation
#define NUM_RUNS (10000 * N)

// the size of the step
#define dt (1.0)

#define TOTAL_TIME (dt * N)

#define SEARCHRANGE (0.7)

#define GRAPHSTEP 10000

// the discretized x values as a function of tau
float x_tau[N];

// returns a random double in range val - step to val + step
float random_val(float val, float step) {
    float normed_rand = ((float)rand()) / ((float)RAND_MAX);
    float scaled_rand = normed_rand * step * 2.0;
    float offset = scaled_rand - step;
    return val + offset;
}

float avg_energy(float * path) {
    float total = 0.0;
    for(size_t idx = 0; idx < N; idx ++) {
        float x_cur = path[idx];
        float pe = 0.5 * m * w * w * x_cur * x_cur;

        float x_plus = path[(idx + 1) % N];
        float dx_plus = x_plus - x_cur;
        float x_minus = path[(idx - 1 + N) % N];
        float dx_minus = x_cur - x_minus;
        float dx_avg = 0.5 * (dx_plus + dx_minus);
        float v_est = dx_avg/dt;
        float ke = 0.5 * m * v_est * v_est;
        total += ke + pe;
    }

    return total/TOTAL_TIME;

}

float analytic_energy(size_t n) {
    return ((float)(2 * n + 1)) * HBAR * w * 0.5;
}

int uniform_distribution(int rangeLow, int rangeHigh) {
    int range = rangeHigh - rangeLow +
                1; //+1 makes it [rangeLow, rangeHigh], inclusive.
    int copies =
        RAND_MAX / range; // we can fit n-copies of [0...range-1] into RAND_MAX
    // Use rejection sampling to avoid distribution errors
    int limit = range * copies;
    int myRand = -1;
    while (myRand < 0 || myRand >= limit) {
        myRand = rand();
    }
    return myRand / copies + rangeLow; // note that this involves the high-bits
}

float tau(uint32_t idx, float temperature) {
    float step = temperature / N;
    return idx * step;
}

size_t accept = 0;
size_t fail = 0;
// runs markov iteration for specific tau coordinate in the chain
float markov_step(float xtau, float temperature, int i) {
    float proposed[N];
    memcpy(proposed, x_tau, N * sizeof(float));
    float xtau_prime = random_val(xtau, SEARCHRANGE);
    assert(!isnan(xtau_prime));
    proposed[i] = xtau_prime;
    float energy_prime = avg_energy(proposed);
    float energy_curr = avg_energy(x_tau);

    float prob_accept = exp( (energy_curr - energy_prime) / (k * temperature));
    if (uniform_distribution(1, 10000) <= prob_accept * 10000) {
        accept += 1;
        return xtau_prime;
    } else {
        fail += 1;
        return xtau;
    }
}

// runs one monte carlo step on the whole polymer chain
void monte_carlo_iteration(float temperature) {
    // choose random index in the chain to move
    int i = uniform_distribution(0, N - 1);

    x_tau[i] = markov_step(x_tau[i], temperature, i);
    assert(!isnan(x_tau[i]));
}

float avg(float *data, size_t len) {
    if (len <= 0) {
        return 0.0;
    }
    float total = 0.0;
    for (size_t idx = 0; idx < len; idx++) {
        total += data[idx];
    }
    return total / ((float)len);
}

float allmin(float *data, size_t len) {
    if (len <= 0) {
        return 0.0;
    }
    float retval = 100000.0;
    for (size_t idx = 0; idx < len; idx++) {
        if (retval > data[idx]) {
            retval = data[idx];
        }
    }
    return retval;
}

void print_state(int j, float * path, FILE * file) {
    for(size_t idx = 0; idx < N; idx ++) {
        float tau = ((float) idx) * dt; 
        float x = path[idx];
        fprintf(file, "%d,%f,%f\n",j, tau, x);
    }
}

int main() {

    float eTotal = 0;
    // number of different simulations
    float *e_chain = (float *)calloc(NUM_RUNS, sizeof(float));
    FILE * ground_state_probability = fopen("ground_state_probability.csv", "w");

    for (int j = 0; j < NUM_RUNS; j++) {
        monte_carlo_iteration(TEMP);
        eTotal = avg_energy(x_tau);
        e_chain[j] = eTotal;
        if(j % GRAPHSTEP == 0) {
            print_state(j, x_tau, ground_state_probability);
        }
    }
    printf("Ground state energy: %f\n", avg(e_chain, NUM_RUNS));
}
