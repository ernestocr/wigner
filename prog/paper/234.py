from sage.all import *

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
    return exp(I * 2 * pi * int(x.trace()) / p)

# \alpha = 0
# P_f = F (Fourier matrix)
F = zero_matrix(SR, d, d)
for i, m in enumerate(FF):
    for j, n in enumerate(FF):
        F[i,j] = omega(m * n) / sqrt(d)
p1 = F

def Pn(diag):
    m = diagonal_matrix(diag)
    return F * m * F.conjugate_transpose()

# \beta = 0
p2 = identity_matrix(d)

# \beta = f(\alpha) = \sigma^6 \alpha + \sigma^3 \alpha^2 + \sigma^5 \alpha^4
p3 = Pn([1, -1, I, I, I, 1, 1, -I])

# \beta = f(\alpha) = \sigma^2 \alpha + \sigma^5 \alpha^2 + \sigma^6 \alpha^4
p4 = Pn([1, -1, -I, 1, I, I, I, 1])

# \beta = f(\alpha) = \sigma^4 \alpha + \sigma^3 \alpha^2 + \sigma^5 \alpha^4
p5 = Pn([1, I, -1, I, -I, I, 1, 1])

# \beta = f(\alpha) = \sigma^3 \alpha
p6 = Pn([1, -I, I, 1, -1, I, 1, I])

# \beta = f(\alpha) = \sigma^5 \alpha + \sigma^5 \alpha^2 + \sigma^6 \alpha^4
p7 = Pn([1, -I, -1, 1, -I, 1, I, -I])

# \beta = f(\alpha) = \sigma \alpha + \sigma^2 \alpha^2 + \sigma \alpha^4
p8 = Pn([1, I, -I, I, 1, 1, I, -1])

# \beta = f(\alpha) = \alpha + \sigma^2 \alpha^2 + \sigma \alpha^4
p9 = Pn([1, 1, 1, I, -1, I, I, -I])

