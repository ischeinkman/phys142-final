#include <float.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include "hermite_polynomial.h"

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

//number of iterations for one markov simulation
#define N (10000)

//lower and upper x-range bounds
#define x0 (-3.0)
#define xD (3.0)

//the number of discrete points in the interval [x0,xD]
#define NUM_POINTS (100)
#define xstep (xD-x0)/NUM_POINTS

//discretized interval of probability distribution
float xvals [NUM_POINTS + 1];

//corresponding probability distribution vals
float prob_distr [NUM_POINTS + 1];

// Lagrangian for the HO 
float L(float x, float v) {
    return 0.5 * m * (v * v - w * w * x * x);
}

// Energy for a given state
// Taken from Wikipedia 
float energy(uint32_t n) {
    return ((float) (2 * n + 1)) * 0.5 * HBAR * w;
}

//returns S(x)
float Sx(uint32_t n, float temperature) {
    return energy(n)/(temperature*k);
}

int uniform_distribution(int rangeLow, int rangeHigh)
{
    int range = rangeHigh - rangeLow + 1; //+1 makes it [rangeLow, rangeHigh], inclusive.
    int copies=RAND_MAX/range; // we can fit n-copies of [0...range-1] into RAND_MAX
    // Use rejection sampling to avoid distribution errors
    int limit=range*copies;    
    int myRand=-1;
    while( myRand<0 || myRand>=limit){
        myRand=rand();   
    }
    return myRand/copies+rangeLow;    // note that this involves the high-bits
}

// The probability to find a particular energy at a 
// particular temperature. Taken from Assignment 3. 
float probability(uint32_t n, float temperature) {
    float power = -1.0*energy(n)/(k * temperature);
    float Z = 1.0; //TODO: Do we even need this?
    return exp(power)/Z; 
}

float hermite(float xvalue, int num) {
    double * as;
    double xv [1] = {xvalue};
    as = h_polynomial_value(1, num, xv);
    return as[0];
}

// N = 0 wavefunction taken from Wikipedia
// C was calced via sagemath
float analytic_state(float x, float n) {
    float fact = 1;
    for (int i = 1; i <= n; i++) {
        fact = fact * i;
    }
    float C = 1.0/sqrt(pow(2.0,n) * fact) * sqrt(sqrt(m * w/(M_PI*HBAR)));
    float power = -0.5 * m * w * x * x/HBAR;
    return C * exp(power) * hermite(x, (int) n);
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

//for generating ground state distribution through markov chain
float prob_trial(float n, float temperature) {
    float n_prime  = (rand() % 3) + (n-1); //generates either n-1, n, or n+1
    float sx_prime = Sx(n_prime, temperature);
    float sx       = Sx(n, temperature);
    if(n_prime < 0)
        return n;
    if(sx_prime < sx) {
        return n_prime;
    }
    else if(sx_prime > sx) {
        float prob_accept = exp(-1.0*sx_prime)/exp(-1.0*sx);

        if(uniform_distribution(1, 10000) <= prob_accept*10000)
            return n_prime;
    }

    return n;
}

void run_quantum_sim(float temperature,
    float (*expect_val)(uint32_t)) {
    float etotal = 0;
    float n = 1; //start from grounded state

    for(int i = 0; i < N; i++) {
        etotal += expect_val(n);
        n = prob_trial(n, temperature);
	for(int j = 0; j < NUM_POINTS + 1; j++) {
            float temp = analytic_state(xvals[j], n);
            prob_distr[j] += temp*temp * probability((uint32_t)n, temperature);
        }
    }
    printf("Ground state energy analytic: %f\n", energy(0));
    printf("Ground state energy approximation: %f\n", etotal/N);
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

void printProbDistr () {
    FILE * file = fopen("prob_distr.txt", "w");
    for(int i = 0; i < NUM_POINTS + 1; i++) {
        char num[32];
	snprintf(num, 32, "%f%12.6f", xvals[i], prob_distr[i]);
	fprintf(file, "%s\n", num);
    }
    fclose(file);
}

int main() {
    //initialize xvals
    for(int i = 0; i < NUM_POINTS + 1; i++) {
        xvals[i] = x0 + i * xstep;
    }
    srand (time(0));
    //choose arbitrary low temperature
    run_quantum_sim(0.1, energy);
    printProbDistr();
    return 0;
}

