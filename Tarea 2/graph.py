import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

df = pd.read_csv('res_erk4.csv')

f1 = plt.figure(1)
plt.plot(df['phi1'], df['phi2'], linewidth=.2, color='black')
plt.xlabel('phi1')
plt.ylabel('phi2')


f2 = plt.figure(2)
plt.plot(df['p1'], df['p2'], linewidth=.2, color='black')
plt.xlabel('p1')
plt.ylabel('p2')


plt.show()
