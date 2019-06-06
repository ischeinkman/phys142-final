#include <float.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

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


int main() {
    return 0;
}

