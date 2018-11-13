import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from math import floor

SPEED = 2

df = pd.read_csv('res_explicit.csv')
fig = plt.figure(figsize=(8, 4), constrained_layout=True)
gs = fig.add_gridspec(2, 2)
fig.suptitle(r'One Dimensional Time-dependent Heat Equation')

ax = fig.add_subplot(gs[:, 0])
ax2 = fig.add_subplot(gs[:, 1])

x = list(map(lambda t: float(t[2:]),list(df)[1:]))

ax2.set_title(r'Explicit Runge-Kutta')
ax2.set_xlabel(r'$x$')
ax2.set_ylabel(r'$T$')
line, = ax2.plot([], [], 'k-', lw=2)
ax2.set_xlim(0, 10)
ax2.set_ylim(-.2, 3)
ax2.grid(True)

time_template = 'steps = %d'
time_text = ax.text(0.5, 0.5, '', transform=ax.transAxes)

def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text

def animate(i):
	i *= SPEED
	temp = df.iloc[i].values.tolist()[1:]
	line.set_data(x, temp)
	time_text.set_text(time_template % df['t'][i])
	return line, time_text

ani = animation.FuncAnimation(fig, animate, np.arange(1, len(df['t']) // SPEED),
                              interval=1, blit=True, init_func=init)

# ani.save('two-body.mp4', fps=24)
plt.show()