import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab
import numpy as np
import pandas as pd 

data=pd.read_csv('ground_state_probability.csv',header=None).values.transpose()
print(len(data))
print(data[0])
print(data[1])
print(data[2])
x_values = data[2]
print(x_values)
num_bins = 2000
n, bins, patches = plt.hist(x_values, num_bins, facecolor = 'blue', alpha = 0.5)
plt.savefig('ground_state_distribution.png', dpi=800)