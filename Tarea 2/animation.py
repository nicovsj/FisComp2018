import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

SPEED = 1

df = pd.read_csv('res_erk4.csv')
fig = plt.figure(figsize=(8, 6), constrained_layout=True)
gs = fig.add_gridspec(2, 3)
fig.suptitle('Péndulo doble')

ax = fig.add_subplot(gs[:, 1:])
ax2 = fig.add_subplot(gs[0, 0])
ax3 = fig.add_subplot(gs[1, 0])

ax2.set_title(r'$\varphi$-space')
ax2.set_xlim(min(df['phi1']), max(df['phi1']))
ax2.set_ylim(min(df['phi2']), max(df['phi2']))
ax2.set_xlabel(r'$\varphi_1$')
ax2.set_ylabel(r'$\varphi_2$')

ax3.set_title(r'$p$-space')
ax3.set_xlim(min(df['p1']), max(df['p1']))
ax3.set_ylim(min(df['p2']), max(df['p2']))
ax3.set_xlabel(r'$p_1$')
ax3.set_ylabel(r'$p_2$')

ax.set_title(r'$(x,z)$-space')
line, = ax.plot([], [], 'ko-', lw=2)
line2, = ax.plot([], [], 'k-', lw=.25)
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$z$')
ax.set_xlim(-2, 2)
ax.set_ylim(-.2, 2.2)
ax.grid(True)

line3, = ax2.plot([], [], 'k-', lw=.25)
line4, = ax3.plot([], [], 'k-', lw=.25)

time_template = 'time = %.1f s'
time_text = ax.text(0.45, 0.2, '', transform=ax.transAxes)

def init():
    line.set_data([], [])
    line2.set_data([], [])
    line3.set_data([], [])
    line4.set_data([], [])
    time_text.set_text('')
    return line, time_text

def animate(i):
    i *= SPEED

    # Seteamos las masas
    thisx = [0, df['x1'][i], df['x2'][i]]
    thisz = [2, df['z1'][i], df['z2'][i]]
    line.set_data(thisx, thisz)

    line2.set_data(df['x2'][:i], df['z2'][:i])

    line3.set_data(df['phi1'][:i], df['phi2'][:i])
    line4.set_data(df['p1'][:i], df['p2'][:i])

    time_text.set_text(time_template % df['t'][i])
    return line, time_text, line2, line3, line4

ani = animation.FuncAnimation(fig, animate, np.arange(1, len(df['x1']) // SPEED),
                              interval=1, blit=True, init_func=init)

# ani.save('two-body.mp4', fps=24)
plt.show()
