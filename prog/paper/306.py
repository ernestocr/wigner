from sage.all import *
import numpy as np

# Three qubits
p = 2
n = 3
d = p**n

# Galois field
f = GF(d, 'x')
# Andres used the default ordering which
# is NOT by powers.
# x = f.gen()
# FF = [f(0)] + [x**j for j in range(d-1)]
FF = f

def omega(x):
    return np.exp(1j * 2 * pi * int(x.trace()) / p)

# P_f = F (Fourier matrix)
# F = zero_matrix(SR, d, d)
F = np.zeros((d,d), dtype='complex128')
for i, m in enumerate(FF):
    for j, n in enumerate(FF):
        F[i,j] = omega(m * n) / np.sqrt(d)
p1 = F # Vertical line

def Pn(diag):
    # m = diagonal_matrix(diag)
    # return F * m * F.conjugate_transpose()
    m = np.diag(diag)
    return F @ m @ F.conj().T

# \beta = 0
p2 = identity_matrix(d) # Horizontal line
# \beta = f(\alpha) = \sigma \alpha
# p3 = Pn([1, -I, I, I, 1, 1, I, -1])
p3 = Pn([1, -1j, 1j, 1j, 1, 1, 1j, -1])
# \beta = f(\alpha) = \sigma^2 \alpha
# p4 = Pn([1, 1, -I, 1, I, I, I, -1])
p4 = Pn([1, 1, -1j, 1, 1j, 1j, 1j, -1])
# \beta = f(\alpha) = \sigma^3 \alpha
# p5 = Pn([1, -I, I, 1, -1, I, 1, I])
p5 = Pn([1, -1j, 1j, 1, -1, 1j, 1, 1j])
# \beta = f(\alpha) = \sigma^4 \alpha
# p6 = Pn([1, I, 1, I, -I, I, 1, -1])
p6 = Pn([1, 1j, 1, 1j, -1j, 1j, 1, -1])
# \beta = f(\alpha) = \sigma^5 \alpha
# p7 = Pn([1, I, -1, 1, -I, 1, I, I])
p7 = Pn([1, 1j, -1, 1, -1j, 1, 1j, 1j])
# \beta = f(\alpha) = \sigma^6 \alpha
# p8 = Pn([1, -1, -I, I, I, 1, 1, I])
p8 = Pn([1, -1, -1j, 1j, 1j, 1, 1, 1j])
# \beta = f(\alpha) = \alpha
# p9 = Pn([1, -1, -1, I, -1, I, I, -I])
p9 = Pn([1, -1, -1, 1j, -1, 1j, 1j, -1j])

# Save to disk
from utils import saveMUBs, checkMUBs

mubs = [p1, p2, p3, p4, p5, p6, p7, p8, p9]
# print('Checking MUBs...')
# checkMUBs(mubs)
print('Saving MUBs...')
saveMUBs(mubs, 'mubs-306.npy')
print('MUBs saved!')