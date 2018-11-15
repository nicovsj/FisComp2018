import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from math import floor

SPEED = 1

df1 = pd.read_csv('res_explicit.csv')
df2 = pd.read_csv('res_implicit.csv')

print(df1.shape, df2.shape)
fig = plt.figure(figsize=(8, 4), constrained_layout=True)
gs = fig.add_gridspec(2, 2)
fig.suptitle(r'One Dimensional Time-dependent Heat Equation')

ax1 = fig.add_subplot(gs[:, 0])
ax2 = fig.add_subplot(gs[:, 1])

x1 = list(map(lambda t: float(t[2:]),list(df1)[1:]))
x2 = list(map(lambda t: float(t[2:]),list(df2)[1:]))

ax1.set_title(r'Explicit Runge-Kutta')
ax1.set_xlabel(r'$x$')
ax1.set_ylabel(r'$T$')
line1, = ax1.plot([], [], 'k-', lw=2)
ax1.set_xlim(0, 10)
ax1.set_ylim(-.2, 3)
ax1.grid(True)


ax2.set_title(r'Implicit Crank-Nicolson')
ax2.set_xlabel(r'$x$')
ax2.set_ylabel(r'$T$')
line2, = ax2.plot([], [], 'k-', lw=2)
ax2.set_xlim(0, 10)
ax2.set_ylim(-.2, 3)
ax2.grid(True)

time_template = 'steps = %d'
time_text = ax1.text(0.5, 0.5, '', transform=ax1.transAxes)

def init():
    line2.set_data([], [])
    line1.set_data([], [])
    time_text.set_text('')
    return line1, line2, time_text

def animate(i):
	i *= SPEED
	temp1 = df1.iloc[i].values.tolist()[1:]
	temp2 = df2.iloc[i].values.tolist()[1:]
	line1.set_data(x1, temp1)
	line2.set_data(x2, temp2)
	time_text.set_text(time_template % df1['t'][i])
	return line1, line2, time_text

ani = animation.FuncAnimation(fig, animate, np.arange(1, len(df2['t']) // SPEED),
                              interval=1, blit=True, init_func=init)

# ani.save('two-body.mp4', fps=24)
plt.show()