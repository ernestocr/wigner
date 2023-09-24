from sage.all import *

"""
This code implements the Wigner function
following Muñoz's paper. It is implemented in two
equivalent ways, using the displacement operators 
and using Wootters' kernel.
"""

N = 3 # number of qubits
F = GF(2**N, 'x') # defaults to Conway's polynomial
x = F.gen()

Id = identity_matrix(SR, 2**N)

# Projection operator (expects a row vector)
def Proj(u, v=None):
    if not v:
        v = u
    u = matrix(u).transpose()
    v = matrix(v).transpose()
    return u.tensor_product(v.conjugate_transpose())

# fix the order of the field
def toInt(k):
    return list(F).index(k)

# Set the additive group character
def chi(k):
    return exp(2 * pi * I * int(k.trace()) / 2)

# Define the Fourier matrix which is
# the eigenbasis of the X_\alpha operators
def Fourier():
    s = zero_matrix(SR, 2**N, 2**N)
    for i, a in enumerate(F):
        for j, b in enumerate(F):
            s[i,j] = chi(a * b) / sqrt(2**N)
    return s
FF = Fourier()

# Define the Z and X operators
def Z(a):
    return diagonal_matrix([chi(a * k) for k in F])

def X(b):
    return FF.conjugate_transpose() * Z(b) * FF

# The following functions depend upon a choice of
# basis, will hard code self-dual basis for ease of use.

N2BASIS = [x,x+1]
N3BASIS = [x**3, x**5, x**6]
# N5BASIS
if N == 2:
    BASIS = N2BASIS
elif N == 3:
    BASIS = N3BASIS

def components(k):
    return vector([(k * m).trace() for m in BASIS])

# Hadamard metric
def h(k):
    return sum([int(bit) for bit in components(k)])

# Muñoz rotation coefficients
def c(alpha, xi, p=1):
    return (-I)**(h(alpha**p * sqrt(xi)**p))

# Phase defined by the rotation coeffs
def phi(tau, nu, p=1):
    t = type(F(0))
    if type(tau) != t:
        tau = F(tau)
    if type(nu) != t:
        nu = F(nu)
        
    if tau == F(0):
        return 1
    return c(tau, tau**-1 * nu, p)

# Displacement operators with phase
def D(a, b, p=1):
    return phi(a, b, p) * Z(a) * X(b)

# Rotation operators
def V(xi, p=1):
    s = zero_matrix(SR, 2**N, 2**N)
    for i, k in enumerate(F):
        s += c(k, xi, p) * Proj(FF[:,i])
    return s

# Displacement kernel
def Delta(a, b):
    s = zero_matrix(SR, 2**N, 2**N)
    for gamma in F:
        for delta in F:
            s += chi(a * delta + b * gamma) * D(gamma,
                                                delta) / 2**N
    return s

# Wootters' kernel
# Hard coding the MUBs
mubs = [FF] + [V(k, p=1) for k in F]
def Wootters(a, b, MUBs=None):
    if not MUBs:
        MUBs = mubs
    op = Proj(FF[:, toInt(a)])
    for xi in F:
        for nu in F:
            d = int(b == xi * a + nu)
            v = MUBs[toInt(xi)+1][:, toInt(nu)]
            op += d * Proj(vector(v))
    return op - Id

# Wigner function
def Wigner(state, a, b, kernel=Delta):
    return (state * kernel(a, b)).trace()

def WignerMatrix(state, kernel=Delta):
    W = zero_matrix(SR, 2**N, 2**N)
    for i, a in enumerate(F):
        for j, b in enumerate(F):
            W[i, j] = real(
                Wigner(state, a, b, kernel)
            ) /2**N
    # return W.matrix_from_rows(range(2**N-1,-1,-1))
    return W
