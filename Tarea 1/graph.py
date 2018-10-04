import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

df = pd.read_csv('res_explicit.txt')
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(211, projection='polar')
# ax.plot(df['phi'], df['rho'])
line, = ax.plot([], [], 'ko', lw=2)
line2, = ax.plot([], [], 'k-', alpha=.5, lw=1.5)
ax.set_rmax(10)
ax.grid(True)

ax2 = fig.add_subplot(212)
ax2.plot(df['rhox'], df['ueffteo'], color='blue', alpha=.3, lw=1.5)
ax2.plot(df['rhox'], df['E'], color='orange')
line3, = ax2.plot([], [], 'bo', lw=2)
line4, = ax2.plot([], [], 'b-', lw=2)
ax2.set_ylabel('Ueff(rho)')
ax2.set_xlabel('rho')
ax2.set_xlim(0, 100)
ax2.set_ylim(min(df['ueffteo'])-1, 40)

time_template = '%.1f years'
time_text = ax.text(0, 0, '', transform=ax.transAxes)


def init():
    line.set_data([], [])
    line2.set_data([], [])
    line3.set_data([], [])
    line4.set_data([], [])
    time_text.set_text('')
    return line, time_text


def animate(i):
    i *= 100
    thisx = [df['phi'][i]]
    thisy = [df['rho'][i]]

    line2.set_data(df['phi'][:i], df['rho'][:i])
    line.set_data(thisx, thisy)
    line3.set_data([df['rho'][i]], [df['ueff'][i]])
    line4.set_data(df['rho'][:i], df['ueff'][:i])
    time_text.set_text(time_template % df['t'][i])
    return line, time_text, line2, line3, line4

ani = animation.FuncAnimation(fig, animate, np.arange(1, len(df['rho'])),
                              interval=1, blit=True, init_func=init)

# ani.save('two-body.mp4', fps=24)
plt.show()
