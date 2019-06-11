import matplotlib
matplotlib.use('Agg')
import pylab
import numpy as np

data=pylab.loadtxt("prob_distr.txt")

pylab.clf()
matplotlib.pyplot.scatter(data[:,0], data[:,1])
pylab.title("Probability Distribution of Ground State of Quantum Harmonic Oscillator")
pylab.xlabel("x")
pylab.ylabel("Probability Density")
pylab.savefig("ground_state_prob_distr.jpg")
pylab.close()

