from sage.all import *
import numpy as np

# Printing options
large_width = 400
np.set_printoptions(linewidth=large_width)

# Three qubits
p = 2
n = 3
d = p**n

# Galois field
f = GF(d, 'x')
FF = f

def omega(x):
    return np.exp(1j * 2 * pi * int(x.trace()) / p)

# \alpha = 0
# P_f = F (Fourier matrix)
# F = zero_matrix(SR, d, d)
F = np.zeros((d, d), dtype='complex128')
for i, m in enumerate(FF):
    for j, n in enumerate(FF):
        F[i,j] = omega(m * n) / np.sqrt(d)

def Pn(diag):
    # m = diagonal_matrix(diag)
    # return F * m * F.conjugate_transpose()
    m = np.diag(diag)
    return F @ m @ F.conj().T

# b = f(a) = x**2 a + x**3 a**2 + x**5 a**4
pf1 = Pn([1,1,1j,1,1j,1j,1j,1])

# b = f(a) = x**6 a**2 + x**3 a**4
pf2 = Pn([1,1,1,1,-1,1,1,-1])

# b = f(a) = a + x**3 a**2 + x**5 a**4
pf3 = Pn([1,-1,1,1j,-1,1j,1j,1j])

# a = g(b) = x**2 b + x**3 b**2 + x**5 b**4
qg1 = np.diag([1,1,1j,1,1j,1j,1j,1])
pg1 = qg1 @ F

# a = g(b) = b + b**6 b**2 + x**3 b**4
qg2 = np.diag([1,-1,-1,1j,1,1j,1j,1j])
pg2 = qg2 @ F

# a = g(b) = x**3 b**2 + x**5 b**4
qg3 = np.diag([1,1,-1,1,1,1,1,-1])
pg3 = qg3 @ F

# a = g(b) = x**6 b + x**3 b**2 + x**5 b**4
qg4 = np.diag([1,-1,1j,1j,1j,1,1,-1j])
pg4 = qg4 @ F

# A = PDP^-1

# Save to disk
# from utils import saveMUBs, checkMUBs

# mubs = [p1, p2, p3, p4, p5, p6, p7, p8, p9]

# print('Saving MUBs...')
# saveMUBs(mubs, 'mubs-162.npy')
# print('MUBs saved!')
# mubs = np.load('mubs-162.npy')
# print('Checking MUBs...')
# checkMUBs(mubs)
