import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np
import pandas as pd
import scipy as sc


with open("ptracks.lst") as fp:
	r = fp.read().split('&')
	hist = r[2].split('\n')[1:-1]
	hist = [[float(x) for x in s.split(',')] for s in hist]

print(hist)


for r in [0.5, 0.6, 0.7, 0.8]:
	Circle((0,0), radius=0.5)

plt.show()



