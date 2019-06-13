#include "hermite_polynomial.h"
#include <assert.h>
#include <float.h>
#include <math.h>
#include <memory.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// Particle mass
#define m (1.0)

// HO Frequency
#define w (1.0)

// Reduced Planck's constant
#define HBAR (1.0)

// Boltzman Constant
#define k (1.0)

// number of discrete points of the polymer chain
#define N (500)

// number of runs of the monte carlo path simulation
#define NUM_RUNS (10000 * N)

#define SEARCHRANGE (0.7)

#define GRAPHSTEP 10000

inline float temp_for_time(float time) { return HBAR / (k * time); }

inline float time_for_temp(float temp) { return HBAR / (k * temp); }

void assert_none_nan(float * path, size_t pathlen) {
    for(size_t idx = 0; idx < pathlen ; idx++) {
        if(isnan(path[idx])) {
            fprintf(stderr, "Got nan at idx: %zu\n", idx);
            assert(!isnan(path[idx]));
        }
    }
}

// returns a random double in range val - step to val + step
float random_val(float val, float step) {
    float normed_rand = ((float)rand()) / ((float)RAND_MAX);
    float scaled_rand = normed_rand * step * 2.0;
    float offset = scaled_rand - step;
    return val + offset;
}


float avg_energy(float *path, size_t pathlen, float total_time) {
    float total = 0.0;
    float dt = total_time / pathlen;
    for (size_t idx = 0; idx < pathlen; idx++) {
        float x_cur = path[idx];
        float pe = 0.5 * m * w * w * x_cur * x_cur;

        float x_plus = path[(idx + 1) % pathlen];
        float dx_plus = x_plus - x_cur;
        float x_minus = path[(idx - 1 + pathlen) % pathlen];
        float dx_minus = x_cur - x_minus;
        float dx_avg = 0.5 * (dx_plus + dx_minus);
        float v_est = dx_avg / dt;
        float ke = 0.5 * m * v_est * v_est;
        total += ke + pe;
    }

    return total / total_time;
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

size_t accept = 0;
size_t fail = 0;

// runs one monte carlo step on the whole polymer chain
void monte_carlo_iteration(float *path, size_t pathlen, float temperature) {
    // choose random index in the chain to move
    int i = uniform_distribution(0, pathlen - 1);
    float xtau = path[i];
    assert(!isnan(xtau));
    float xtau_prime = random_val(xtau, SEARCHRANGE);
    assert(!isnan(xtau_prime));
    float cur_energy = avg_energy(path, pathlen, time_for_temp(temperature));
    path[i] = xtau_prime;
    float prime_energy = avg_energy(path, pathlen, time_for_temp(temperature));
    float prob_accept = exp((cur_energy - prime_energy) / (k * temperature));
    int scale_factor = 10000;
    if (uniform_distribution(1, scale_factor) <=
        ((int)scale_factor * prob_accept)) {
        accept += 1;
        path[i] = xtau_prime;
    } else {
        fail += 1;
        path[i] = xtau;
    }
    assert(!isnan(path[i]));
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

void print_state(int j, float dt, float *path, size_t pathlen, FILE *file) {
    for (size_t idx = 0; idx < pathlen; idx++) {
        float tau = ((float)idx) * dt;
        float x = path[idx];
        fprintf(file, "%d,%f,%f\n", j, tau, x);
    }
}

float simulate_temperature(float *path, size_t pathlen, float temperature,
                           int runs) {
    float eTotal = 0;
    float *e_chain = (float *)calloc(runs, sizeof(float));

    for (int j = 0; j < runs; j++) {
        monte_carlo_iteration(path, pathlen, temperature);
        eTotal = avg_energy(path, pathlen, time_for_temp(temperature));
        e_chain[j] = eTotal;
    }
    float retval = avg(e_chain, runs);
    free(e_chain);
    return retval;
}

#define TESTEP 0.001
#define TESTART (1.0 / 500.0)
#define TESTOP 5.0
#define TESIZE ((TESTOP - TESTART + 1.0) / TESTEP)
#define TERUNS 10000

void energy_plot() {
    FILE *output = fopen("energy_for_temp.csv", "w");
    for (float temp = TESTART; temp < TESTOP; temp += TESTEP) {
        float *path = (float *)calloc(N, sizeof(float));
        float cur_e = simulate_temperature(path, N, temp, TERUNS);
        fprintf(output, "%f,%f\n", temp, cur_e);
        free(path);
    }
    fclose(output);
}

int main() {

    const float ground_state_temp = 1.0/500.0;
    const float ground_state_time = time_for_temp(ground_state_temp);
    const float ground_state_dt = ground_state_time/((float) N);

    float *e_chain = (float *)calloc(NUM_RUNS, sizeof(float));

    float *x0xtau = (float *)calloc(NUM_RUNS, sizeof(float));
    float *x0xtaup1 = (float *)calloc(NUM_RUNS, sizeof(float));

    float *path = (float *) calloc(N, sizeof(float));

    FILE *ground_state_probability = fopen("ground_state_probability.csv", "w");

    for (int j = 0; j < NUM_RUNS; j++) {
        assert_none_nan(path, N);
        monte_carlo_iteration(path, N, ground_state_temp);
        e_chain[j] = avg_energy(path, N, time_for_temp(ground_state_temp));
        if (j % GRAPHSTEP == 0) {
            print_state(j / GRAPHSTEP, ground_state_dt, path, N, ground_state_probability);
        }
        x0xtaup1[j] = path[0] * path[255];
        x0xtau[j] = path[0] * path[245];
    }
    float top = avg(x0xtaup1, NUM_RUNS);
    float bottom = avg(x0xtau, NUM_RUNS);
    float dt = time_for_temp(ground_state_temp) / N;
    float diff = -1.0 / (dt * 10.0) * log(top / bottom) / log(exp(1));
    float E1 = analytic_energy(0) + diff;
    printf("Ground state energy: %f\n", avg(e_chain, NUM_RUNS));
    printf("First excited state energy: %f\n", E1);
    fclose(ground_state_probability);
    energy_plot();
}