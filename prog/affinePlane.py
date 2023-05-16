# %%
from sage.all import *
import numpy as np
import matplotlib.pyplot as plt
# %%
def affineAx(F):
    d = F.order
    fig, ax = plt.subplots()

    ax.xaxis.set_ticks_position('bottom')
    # ax.set_xticks(range(d))
    # ax.set_yticks(range(d))

    # with F.repr('power'):
    #     print(F.elements)
    #     labels = ['$' + str(e) + '$' for e in F]
    #     ax.set_xticklabels(labels)
    #     ax.set_yticklabels(labels)
    
    return fig, ax
# %%
def DesarguesSpread(F):
    lines = [np.array([(0,u) for u in F.elements])]
    for m in F.elements:
        lines.append(
            np.array([(u, m*u) for u in F.elements]))
    return lines
# %%
F = GF(9, repr='int')
# %%
fig, ax = affineAx(F)

lines = DesarguesSpread(F)

for i, line in enumerate(lines):
    alpha = 0.2
    if i in [0,1,2,4]:
        alpha = 0.8

    ax.plot(
        line[:,0], line[:,1],
        marker='o',
        linewidth=0.2,
        markersize=6,
        alpha=alpha
    )
# %%
