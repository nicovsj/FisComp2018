import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from math import floor, ceil

df = pd.read_csv('res_explicit.csv')

x = list(map(lambda t: float(t[2:]),list(df)[1:]))
temp = df.iloc[30].values.tolist()[1:]

print(x)

print(floor(1/0.5))

plt.plot(x, temp)
plt.show()