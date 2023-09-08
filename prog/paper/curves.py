#%%
from sage.all import *
import matplotlib.pyplot as plt

#%%
# Let's plot the curves!
p = 2
n = 3
d = p**n

F = GF(d, 'x')
x = F.gen()

#%%
def curveIntegers(curve):
    return [(toInt(p[0]), toInt(p[1])) for p in curve]

def plotCurves(curves, title=None):
    fig, ax = plt.subplots()
    labels = ['oo'] + [str(k) for k in F]
    for i, curve in enumerate(curves):
        ax.plot(*zip(*curveIntegers(curve)), marker='o', label=i) 
    plt.title(title)
    plt.legend(loc='upper right')
    plt.show()

def plotCurve(curve):
    fig, ax = plt.subplots()
    ax.plot(*zip(*curveIntegers(curve)), marker='o')
    ax.grid()
    plt.show()
# %%
curve = [
    (F(0), F(0)),
    (F(0), x**5),
    (F(1),x),
    (F(1),x**6),
    (x,x**2),
    (x,x**3),
    (x**2,x**2),
    (x**2,x**3),
    (x**3,x**4),
    (x**3,F(1)),
    (x**4,F(0)),
    (x**4,x**5),
    (x**5,x),
    (x**5,x**6),
    (x**6,x**4),
    (x**6,F(1))
]
plotCurve(curve)
iso = {}
for c in curve:
    for b in curve:
        if (c[0]*b[1] - c[1]*b[0]).trace() == 0:
            iso.add(c)
print(iso)
# %%
# (3,0,6) (Standard)

curves = []

# \alpha = 0
vertical = [(0, k) for k in F]
curves.append(vertical)

# \beta = 0
horizontal = [(k, 0) for k in F]
curves.append(horizontal)

def toInt(k):
    return list(F).index(k)

# \beta = \sigma \alpha
for j in range(1,7):
    curves.append([(k, (x**(j)) * k) for k in F])
curves.append([(k, k) for k in F])

plotCurves(curves)

# %%
# (1,6,2)
curves = []

# \alpha = 0
vertical = [(0, k) for k in F]
curves.append(vertical) 

# \beta = \sigma \alpha^2 + \sigma \alpha^4
curves.append([(k, x * k**2 + x * k**4) for k in F])
# \beta = \sigma \alpha + \sigma \alpha^2 + \sigma \alpha^4
curves.append([(k, x * k + x * k**2 + x * k**4) for k in F])
# \beta = \sigma^2 \alpha + \sigma \alpha^2 + \sigma \alpha^4
curves.append([(k, x**2 * k + x * k**2 + x * k**4) for k in F])
# \beta = \sigma^3 \alpha + \sigma \alpha^2 + \sigma \alpha^4
curves.append([(k, x**3 * k + x * k**2 + x * k**4) for k in F])
# \beta = \sigma^4 \alpha + \sigma \alpha^2 + \sigma \alpha^4
curves.append([(k, x**4 * k + x * k**2 + x * k**4) for k in F])
# \beta = \sigma^5 \alpha + \sigma \alpha^2 + \sigma \alpha^4
curves.append([(k, x**5 * k + x * k**2 + x * k**4) for k in F])
# \beta = \sigma^6 \alpha + \sigma \alpha^2 + \sigma \alpha^4
curves.append([(k, x**6 * k + x * k**2 + x * k**4) for k in F])
# \beta = \alpha + \sigma \apha^2 + \sigma \alpha^4
curves.append([(k, k + x * k**2 + x * k**4) for k in F])

plotCurves(curves)
# %%
# (2,3,4)

curves = []

# \beta = 0
horizontal = [(k, 0) for k in F]
curves.append(horizontal)

# \alpha = 0
vertical = [(0, k) for k in F]
curves.append(vertical) 

# \beta = \sigma^6 \alpha + \sigma^3 \alpha^2 + \sigma^5 \alpha^4
curves.append([(k, x**6 * k + x**3 * k**2 + x**5 * k**4) for k in F])
# \beta = \sigma^2 \alpha + \sigma^5 \alpha^2 + \sigma^6 \alpha^4
curves.append([(k, x**2 * k + x**5 * k**2 + x**6 * k**4) for k in F])
# \beta = \sigma^4 \alpha + \sigma^3 \alpha^2 + \sigma^5 \alpha^4
curves.append([(k, x**4 * k + x**3 * k**2 + x**5 * k**4) for k in F])
# \beta = \sigma^3 \alpha
curves.append([(k, x**3 * k) for k in F])
# \beta = \sigma^5 \alpha + \sigma^5 \alpha^2 + \sigma^6 \alpha^4
curves.append([(k, x**5 * k + x**5 * k**2 + x**6 * k**4) for k in F])
# \beta = \sigma \alpha + \sigma^2 \alpha^2 + \sigma \alpha^4
curves.append([(k, x * k + x**2 * k**2 + x * k**4) for k in F])
# \beta = \alpha + \sigma^2 \alpha^2 + \sigma \alpha^4
curves.append([(k, k + x**2 * k**2 + x * k**4) for k in F])
# %%
plotCurves(curves, title='(3,4,2)')
# %%
for curve in curves:
    for p1 in curve:
        for p2 in curve:
            p = (p1[0] + p2[0], p1[1] + p2[1])
            if p not in curve:
                raise Exception('Not closed under addition!')
# %%
