# %%
import numpy as np
import matplotlib.pyplot as plt
import galois

d = 2
n = 3
GF = galois.GF(d**n, irreducible_poly='x^3 + x^2 + 1')
GF.repr('power')

# defaults to integer ordering, lets change to power order
# is this important?
omega = GF.primitive_element
F = omega**np.arange(0, GF.order - 1)
F = np.insert(F, 0, 0)
# F = GF.elements
# %%
def heat(m, ax, labels=True):
    m = np.rot90(m, axes=(1,0)).T
    ax.matshow(
        m,
        cmap='viridis',
        origin='lower'
    )
    if labels:
        ax.xaxis.set_ticks_position('bottom')
        ax.set_xticks(range(d**n))
        ax.set_yticks(range(d**n))
        ax.set_xticklabels(F)
        ax.set_yticklabels(F)
    else:
        ax.axis('off')
    return ax
# %%
def embedLine(line):
    N = d**n
    grid = np.zeros((N,N))
    for p in np.array(line):
        grid[p[0], p[1]] = 1.0
    return np.rot90(grid)

fig, ax = plt.subplots()
heat(embedLine(GF([(k,F[2]*k) for k in F])), ax)
# %%
fig, axes = plt.subplots(d**n + 1, d**n, sharey=True)
for (j,c) in enumerate(F):
    vertical_line = GF([(c,k) for k in F])
    heat(embedLine(vertical_line), axes[0][j], labels=False)

    for (i,m) in enumerate(F):
        line = GF([(k, m * k + c) for k in F])
        heat(embedLine(line), axes[i+1][j], labels=False)
# %%
