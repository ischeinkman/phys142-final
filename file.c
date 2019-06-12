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

//number of discrete points of the polymer chain
#define N (500)

//number of runs of the monte carlo path simulation
#define NUM_RUNS (1000)

//the size of the step
#define delta (2)

//the discretized tau values of the polymer chain
float tau_vals[N];

//the discretized x values as a function of tau
float x_tau[N];

//returns a random double in range val - step to val + ste[
float random_val (float val, float step) {
    return (step + step) * ( (float) rand() / (float) RAND_MAX) + val - step;
}

// Energy for a given state
// Taken from Wikipedia
float chain_energy(float xtau, int i) {
    int index_plus = (i+1)%N;
    int index_minus = (i-1)%N;
    float xi_plus = x_tau[index_plus] - xtau;
    float xi_minus = x_tau[index_minus] - xtau;
    return delta * (0.5*m*( (xi_plus*xi_plus) + (xi_minus*xi_minus) ) / (delta * delta) +
        0.5 * m * w * w * xtau * xtau);
}

int uniform_distribution(int rangeLow, int rangeHigh) {
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

//initializes tau_array to correct discrete values for given temp
void initialize_values (float temperature) {
    float step = temperature/N;
    for(int i = 0; i <= N; i++) {
        tau_vals[i] = i*step;
    }
}

//runs markov iteration for specific tau coordinate in the chain
float markov_step(float xtau, float temperature, int i) {
    float xtau_prime = random_val(xtau, delta);
    float energy_prime = chain_energy(xtau_prime, i);
    float energy_curr = chain_energy(xtau, i);

    if(energy_curr < energy_prime) {
        return xtau_prime;
    }

    else if (energy_curr > energy_prime) {
        float prob_accept = exp(-energy_curr/temperature)/exp(-energy_prime/temperature);
	if (uniform_distribution(1,10000) <= prob_accept*10000)
            return xtau_prime;
    }
    return xtau;
}

//runs one monte carlo step on the whole polymer chain
void monte_carlo_iteration(float temperature) {
    //choose random index in the chain to move
    int i = (int)((rand() * 1.0) / RAND_MAX * N);

    x_tau[i] = markov_step(tau_vals[i], temperature, i);
}

int main(){
    float T = 1000;
    initialize_values(T);

    float eTotal = 0;
    //number of different simulations
        for(int j = 0; j <2000*N; j++) {
            for(int h = 0; h < N; h++) {
		 float diff = x_tau[(h+1)%N] - x_tau[h];
		 eTotal += T * delta * (0.5*m*(diff * diff) / (delta*delta) +
	             0.5 * m * w * w  * x_tau[h]*x_tau[h]);
            }
            monte_carlo_iteration(T);
	}
    printf("%f\n", eTotal/2000/N);
}
