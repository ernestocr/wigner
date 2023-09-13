import numpy as np
from sage.all import *
from wigner import Wigner, PlotWignerFunction
import matplotlib.pyplot as plt

# ---------------------------------------------------
# Finding stabilizer states

# Three qubits
p = 2
n = 3
d = p**n

# Galois field
FF = GF(d, 'x')
x = FF.gen()

def omega(x):
    return np.exp(1j * 2 * np.pi * int(x.trace()) / p)

F = np.zeros((d, d), dtype='complex128')
for i, m in enumerate(FF):
    for j, n in enumerate(FF):
        F[i,j] = omega(m * n) / np.sqrt(d)

def Z(a):
    return np.diag([omega(a * n) for n in FF])

def X(b):
    return F @ Z(b) @ F.conj().T

def D(a,b):
    return Z(a) @ X(b)

def FindStabilizer(mubs, curves):
    stabilizers = []
    for i in range(8):
        disps = [D(t, curves[i](t)) for t in FF]

        # Commutativity sanity check
        # for disp1, disp2 in zip(disps, disps):
        #     if not np.all(disp1 @ disp2 == disp2 @ disp1):
        #         raise Exception('Operators not commuting.')

        # Numerically compute the eigenvalues
        for disp in disps:
            print(np.round(np.linalg.eigvals(disp)))

        # for j in range(8):
        #     v = mubs[(i+1)*8:(i+2)*8, j]
        #     a = True
        #     for disp in disps:
        #         if not np.all(np.isclose(disp @ v, v)):
        #             print(np.isclose(v, disp @ v))
        #             print(v, disp @ v)
        #             a = False
        #             break
        #     if a:
        #         stabilizers.append((i,j))
    return stabilizers

mubs = np.load('mubs-306.npy')
curves = [
    lambda t: 0,
    lambda t: x * t,
    lambda t: x**2 * t,
    lambda t: x**3 * t,
    lambda t: x**4 * t,
    lambda t: x**5 * t,
    lambda t: x**6 * t,
    lambda t: t
]
w = Wigner(8, mubs)
w.LoadCurves(curves)

# Load Wigner function for 234 mubs
mubs234 = np.load('mubs-234.npy')
# w = Wigner(8, mubs234)
curves = [
    lambda t: 0,
    lambda t: x**6 * t + x**3 * t**2 + x**5 * t**4,
    lambda t: x**2 * t + x**5 * t**2 + x**6 * t**4,
    lambda t: x**4 * t + x**3 * t**2 + x**5 * t**4,
    lambda t: x**3 * t,
    lambda t: x**5 * t + x**5 * t**2 + x**6 * t**4,
    lambda t: x**1 * t + x**2 * t**2 + x**1 * t**4,
    lambda t: t + x**2 * t**2 + x * t**4
]
# results = FindStabilizer(mubs, curves)
# print(len(results))

# # Plot each stabilizer state using the 234 Wigner function
# for i, j in result:
#     w.LoadState(w.Proj(mubs[i*8:(i+1)*8, j]))
#     w.WignerMatrix(recalc=True)
#     fig, ax = PlotWignerFunction(w)
#     plt.show()

stabilizers = []
for i in range(9):
    for j in range(8):
        v = w.Proj(mubs234[8*i:8*(i+1), j])
        w.LoadState(v)
        w.WignerMatrix(recalc=True)
        if np.all(w.cache >= -1e-16):
            stabilizers.append((i,j))