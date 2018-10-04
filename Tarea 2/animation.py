import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

SPEED = 10

df = pd.read_csv('res_erk4.csv')
fig = plt.figure(figsize=(12, 8), constrained_layout=True)
gs = fig.add_gridspec(2, 3)
fig.suptitle('Péndulo doble')

ax = fig.add_subplot(gs[:, 1:])
ax2 = fig.add_subplot(gs[0, 0])
ax3 = fig.add_subplot(gs[1, 0])

ax2.set_xlim(-1.5, 1.5)
ax2.set_ylim(2, -2)
ax2.set_xlabel('phi1')
ax2.set_ylabel('phi2')

ax3.set_xlim(-4.5, 4.5)
ax3.set_ylim(3, -3)
ax3.set_xlabel('p1')
ax3.set_ylabel('p2')

# ax.plot(df['phi'], df['rho'])
line, = ax.plot([], [], 'ko-', lw=2)
line2, = ax.plot([], [], 'k-', alpha=.5, lw=1.5)
ax.set_xlabel('x')
ax.set_ylabel('z')
ax.set_xlim(-2, 2)
ax.set_ylim(-.2, 2.2)
ax.grid(True)

line3, = ax2.plot([], [], 'k-')
line4, = ax3.plot([], [], 'k-')

# ax2 = fig.add_subplot(212)
# ax2.plot(df['rhox'], df['ueffteo'], color='blue', alpha=.3, lw=1.5)
# ax2.plot(df['rhox'], df['E'], color='orange')
# line3, = ax2.plot([], [], 'bo', lw=2)
# line4, = ax2.plot([], [], 'b-', lw=2)
# ax2.set_ylabel('Ueff(rho)')
# ax2.set_xlabel('rho')
# ax2.set_xlim(0, 100)
# ax2.set_ylim(min(df['ueffteo'])-1, 40)

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
