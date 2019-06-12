import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv('energy_for_temp.csv', header=None).values.transpose()
t_data = data[0]
e_data = data[1]

plt.scatter(t_data, e_data)
plt.xlabel('Temperature')
plt.ylabel('Average Energy')
plt.savefig('t_vs_e.png', dpi=500)