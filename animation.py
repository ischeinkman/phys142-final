import matplotlib
matplotlib.use('Agg')
import pandas as pd
import pylab
import numpy as np
import math

num=500
data = pd.read_csv("ground_state_probability.csv").values.transpose()
tauvals = data[1]
xvals = data[2]
i=0
j=0
while (i<len(xvals)):
    pylab.clf()
    pylab.plot(tauvals[i:(i+num - 1)], xvals[i:(i+num - 1)], ms=0.2)

    pylab.title("Markov Chain Path Sampling")
    pylab.xlabel("T")
    pylab.axis([0,500,-4,4])
    pylab.ylabel("x(T)")
    pylab.rcParams['figure.figsize'] = 10, 5
    pylab.savefig("out/" + ('%03d' % j)+ '.png', dpi=800)
    pylab.close()
    j+= 1
    i+= num
