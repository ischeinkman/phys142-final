#include <float.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <stdio.h>

#define min(a, b) ((a > b) ? b : a)

// Particle mass
#define m (1.0)

// HO Frequency
#define w (1.0)

// Reduced Planck's constant
#define HBAR (1.0)

// Boltzman Constant
#define k (1.0)

// The temperature to use for the 
// ground state energy calcs
#define LOW_TEMP (FLT_EPSILON)

// Lagrangian for the HO 
float L(float x, float v) {
    return 0.5 * m * (v * v - w * w * x * x);
}

// Energy for a given state
// Taken from Wikipedia 
float energy(uint32_t n) {
    return ((float) (2 * n + 1)) * 0.5 * HBAR * w;
}

// The probability to find a particular energy at a 
// particular temperature. Taken from Assignment 3. 
float probability(float energy, float temperature) {
    float power = -energy/(k * temperature);
    float Z = 1.0; //TODO: Do we even need this?
    return exp(power)/Z; 
}

// N = 0 wavefunction taken from Wikipedia
// C was calced via sagemath
float analytical_ground(float x) {
    float C = sqrt(2.0 * M_PI * HBAR/(m * w));
    float power = -0.5 * m * w * x * x/HBAR;
    return C * exp(power);
}

float normed_rand() {
    return ((float) rand())/((float) RAND_MAX);
}


//Problem 1: Monte-Carlo chain a point particle in a LOW_TEMP heat bath and HO potential 
#define p1_steps 10000000
#define p1_radius 2.0
#define dt 1.0

float next_energy(float cur_x, float next_x, float time_between) {
    float x_est = (next_x + cur_x)/2.0;
    float pe = 0.5 * m * w * w * x_est * x_est;

    float v_est = (next_x - cur_x)/time_between; 
    float ke = 0.5 * m * v_est * v_est;
    return ke + pe; 
}

float p1() {
    float cur_x = 0.0;
    float cur_E = 0.0;
    float t_here = 0.0;
    for(size_t idx = 0; idx < p1_steps; idx++) {
        t_here += 1.0;
        float next_x = p1_radius * (normed_rand() -0.5) + cur_x;
        float next_E = next_energy(cur_x, next_x, t_here);

        float metropolis_factor = min(1.0, probability(next_E, LOW_TEMP)/probability(cur_E, LOW_TEMP));
        float flag = normed_rand();
        if(flag < metropolis_factor) {
            cur_x = next_x;
            cur_E = next_E;
            t_here = 0.0;
        }
    }
    return cur_E;
}



int main() {
    fprintf(stdout, "Ground state energy: %f\n", p1());
    return 0;
}

