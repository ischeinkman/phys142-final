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
    pylab.plot(tauvals[i:(i+num)], xvals[i:(i+num)], ms = .2)

    pylab.title("Markov Chain Path Sampling")
    pylab.xlabel("T")
    pylab.ylabel("x(T)")
    pylab.savefig("out/out" + ('%03d' % j)+ '.jpg')
    pylab.close()
    j+= 1
    i+= num
