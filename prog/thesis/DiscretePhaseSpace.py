# %%
import numpy as np
import matplotlib.pyplot as plt
from sage.all import *

p = 2
n = 3
d = p**n
F = GF(d, 'x')
x = F.gen()

F = [F(0)] + [x**i for i in range(d-1)]
# %%
def heat(m, ax, labels=True):
    ax.matshow(
        m,
        cmap='viridis',
        origin='lower'
    )
    if labels:
        ax.xaxis.set_ticks_position('bottom')
        ax.set_xticks(range(d))
        ax.set_yticks(range(d))
        ax.set_xticklabels(F)
        ax.set_yticklabels(F)
    else:
        ax.axis('off')
    return ax
# %%
def toInt(x):
    return F.index(x)

def embedLine(line):
    grid = np.zeros((d, d))
    for p in line:
        i, j = toInt(p[0]), toInt(p[1])
        grid[i, j] = 1.0
    return grid.T

fig, ax = plt.subplots()
heat(embedLine([(k,x*x*k) for k in F]), ax)
# %%
fig, axes = plt.subplots(d + 1, d,
                         sharey=True,
                         figsize=(10,10))
for (j,c) in enumerate(F):
    vertical_line = [(c,k) for k in F]
    heat(embedLine(vertical_line), axes[0][j], labels=False)

    for (i,m) in enumerate(F):
        line = [(k, m * k + c) for k in F]
        heat(embedLine(line), axes[i+1][j], labels=False)
# %%
