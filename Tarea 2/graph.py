import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from random import randint


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
k = randint(0, len(df['x2']))
ax.plot([0, df['x1'][k], df['x2'][k]], [2, df['z1'][k], df['z2'][k]],
        'ko-', lw=1)
for _ in range(2):
    k += len(df['x1']) // 300
    ax.plot([0, df['x1'][k], df['x2'][k]], [2, df['z1'][k], df['z2'][k]],
            'ko--', lw=1, alpha=.5)

ax.plot(df['x2'], df['z2'], 'k-', lw=.25)
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$z$')
ax.set_xlim(-2, 2)
m = 2.2
if max(df['z2']) > 2:
    m = 3
ax.set_ylim(-.2, m)
ax.grid(True)

ax2.plot(df['phi1'], df['phi2'], 'k-', lw=.25)
ax3.plot(df['p1'], df['p2'], 'k-', lw=.25)

#plt.savefig('pend.png')
plt.show()
