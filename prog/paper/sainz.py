from sage.all import *
import numpy as np

# Build the field, the Galois ring and the Teichmüller set

p = 2
N = 3 # number of qubits

F = GF(p**N, 'x')
x = F.gen()

R = PolynomialRing(Integers(4), 't')
t = R.gen()

GR = R.quotient(t**3 + 2*t**2 + t - 1, 'w')
w = GR.gen()

# Get field element index
def toInt(k):
    return list(F).index(k)

# Teichmüller set
T = [GR(0)] + [w**j for j in range(1, 2**N - 1)] + [GR(1)]

# Lift a field element to its Teichmüller representative
def TeichLift(k):
    return T[toInt(k)]
hat = TeichLift # shortcut

# Now build the group character and the generalized
# Pauli operators

def chi(k):
    return np.exp(np.pi * 1j * int(k.trace()))

def Proj(u, v=None):
    if v is None:
        v = u
    return np.outer(u, v.conj().T)

Id = np.eye(2**N, 2**N)

def Fourier():
    s = np.zeros((2**N, 2**N), dtype='complex128')
    for i, a in enumerate(F):
        for j, b in enumerate(F):
            s[i,j] = chi(a * b) / np.sqrt(2**N)
    return s
FF = Fourier()

def Z(a):
    return np.diag([chi(a * k) for k in F])

def X(b):
    return FF.conj().T @ Z(b) @ FF

# Sainz's phase, it uses the Galois ring trace
# down to Z_4. SageMath handles the 2-adic
# decomposition to define the trace.
def phi(a, b):
    if not isinstance(a, type(w)):
        a = hat(a)
    if not isinstance(b, type(w)):
        b = hat(b)

    return (1j)**(int((a * b).trace()))

def D(a, b):
    return phi(a, b) * (Z(a) @ X(b))

# Definiion of the rotation operators
# for the Desarguesian spread.
# This function gives the full basis
# corresponding to a given parameter \mu.
def V(mu):
    s = np.zeros((2**N, 2**N), dtype='complex128')
    for i, k in enumerate(F):
        s += phi(k, mu * k) * Proj(FF[:, i])
    return s

# We can obtain the whole set of MUBs, ordered
# by the field elements in the following way:
mubs306 = [FF] + [V(mu) for mu in F]

# Now for the Wigner function.

# Wootters kernel for Desarguesian spread
def Wootters(a, b):
    op = Proj(mubs306[0][:, toInt(a)]) # Fourier basis
    for xi in F:
        for nu in F:
            d = float(b == xi * a + nu)
            op += d * Proj(mubs306[toInt(xi)+1][:,toInt(nu)])
        # nu = b - xi * a
        # op += Proj(mubs306[toInt(xi)+1][:,toInt(nu)])
    return op - Id

# Definition of the Wigner function
def Wigner(state, a, b, kernel):
    return (state @ kernel(a, b)).trace()

def WignerMatrix(state, kernel):
    W = np.zeros((2**N, 2**N))
    for i, a in enumerate(F):
        for j, b in enumerate(F):
            W[i, j] = np.real(Wigner(state, a, b, kernel)) / (2**N)
    return W

def WW(state, approx=True):
    w = WignerMatrix(Proj(state), Wootters)
    if approx:
        return np.round(w, 3)
    return w

def Prob(w, curve):
    s = 0
    for i, a in enumerate(F):
        for j, b in enumerate(F):
            s += w[i,j] * float(b == curve(a))
    return s


# Part of Sainz's proof

mubs162 = np.load('sainz/162.npy')
mubs162 = [mubs162[j*8:(j+1)*8] for j in range(9)]

# alpha = x
# beta  = x
# try:
#     for s in F:
#         curve = lambda t: s * t + t**2 + t**4
#         curve_hat = lambda t: hat(s * t) + hat(t**2) + hat(t**4)
        
#         for k in F:
#             left_side = 0 + 0j
            
#             for mu in F:
#                 nu = beta - mu * alpha

#                 psi_mu_nu = mubs306[toInt(mu)+1][:,toInt(nu)]
#                 psi_f_k = mubs162[toInt(s)+1][:,toInt(k)]

#                 inp = abs(psi_mu_nu.conj() @ psi_f_k)**2
#                 left_side += inp

#                 ss = 0
#                 for xi in F:
#                     c_xif = phi(xi, curve_hat(xi))
#                     c_ximu = phi(xi, xi * mu)
#                     char = chi(xi * (-k + nu))
#                     d = int(mu * xi == curve(xi))
#                     ss += c_xif * conjugate(c_ximu) * char * d / 2**N
                
#                 if not np.isclose(inp, ss):
#                     raise Exception('No equality!', mu, nu, s, k)
                
#             right_side = 1 - 1/(2**N) + sum(
#                 [chi(xi * (beta-curve(alpha)-k)) for xi in F]
#             ) / (2**N)
                
#             print(
#                 np.isclose(
#                     left_side, right_side
#                 ),
#                 np.round(left_side, 5), np.round(right_side, 5),
#                 s, k
#             )
# finally:
#     print('Equation (40) is valid!')


# for mu in F:
#     curve = lambda t: mu * t + t**2 + t**4
#     for alpha in F:
#         for xi in F:
#             if chi(alpha * curve(xi)) != chi(xi * curve(alpha)):
#                 raise Exception('Not a symplectic presemifield.')
# print('Curves satisfy the needed property!')

alpha = x
beta  = x
for s in F:
    curve = lambda t: s * t + t**2 + t**4
    curve_hat = lambda t: hat(s * t) + hat(t**2) + hat(t**4)
    for k in F:
        s40 = 0
        s42 = 1
        for mu in F:
            nu = beta - mu * alpha
            for xi in F:
                c_xif = phi(hat(xi), curve_hat(xi))
                c_ximu = phi(xi, mu * xi)
                char = chi(xi * (-k + beta - mu * alpha))
                d = int(mu * xi == curve(xi))
                s40 += c_xif * conjugate(c_ximu) * char * d / (2**N)
        
        for xi in list(F)[1:]:
            s42 += chi(xi * (-k + beta - curve(alpha))) / (2**N)
        
        print(np.round(s40, 4), np.round(s42, 4))


# TESTING

# Testing the rotation operators

from utils import checkPhasePointOperators

class TestSuite:

    def __init__(self):
        pass

    def run(self):
        self.RotationOperators()
        self.RotationCoefficients()
        self.KernelTest()
        self.TomographicProperties()
        self.Covariance()

        print('All tests have passed successfully!')

    def RotationOperators(self):
        try:
            for mu in F:
                for alpha in F:
                    v_op = V(mu)
                    z_op = Z(alpha)
                    rot_op = v_op @ z_op @ v_op.conj().T
                    disp = D(alpha, mu * alpha)
                    if not np.all(np.isclose(rot_op, disp)):
                        raise Exception('Rotation failed!') 
        finally:
            print('Rotation operators work!')

    def RotationCoefficients(self):
        try:
            for mu in F:
                for k in F:
                    for a in F:
                        c_k = phi(k, k * mu)
                        c_a = phi(a, a * mu)
                        c_ak = phi(k + a, (k + a) * mu)
                        char = chi(mu * a * k)
                        if not np.isclose(c_k * c_a, c_ak * char):
                            raise Exception(
                                'Recurrence relation does not hold!',
                            mu, k, a)
        finally:
            print('Recurrence relation holds for chosen phase!')

    def KernelTest(self):
        ops = []
        for a in F:
            for b in F:
                ops.append(Wootters(a, b))
        checkPhasePointOperators(ops)
        print('Normalized', np.all(np.isclose(sum(ops) / 8, Id)))

    def TomographicProperties(self):
        try:
            # Vertical and horizontal lines
            for nu in F:
                op = np.zeros((2**N, 2**N), dtype='complex128')
                for a in F:
                    op += Wootters(a, nu)
                if not np.all(np.isclose(op / 2**N, Proj(Id[:, toInt(nu)]))):
                    raise Exception('Error with Z eigenbasis!')
                
                op = np.zeros((2**N, 2**N), dtype='complex128')
                for b in F:
                    op += Wootters(nu, b)
                if not np.all(np.isclose(op / 2**N, Proj(FF[:, toInt(nu)]))):
                    raise Exception('Error with X eigenbasis!')
                
            # Arbitrary curve
            for mu in F:
                for nu in F:
                    op = np.zeros((2**N, 2**N), dtype='complex128')
                    for a in F:
                        for b in F:
                            op += Wootters(a, b) * float(b == mu * a + nu)
                    v = mubs306[toInt(mu)+1][:,toInt(nu)]
                    if not np.all(np.isclose(op / 2**N, Proj(v))):
                        raise Exception('Error with Z_a X_f(a) eigenbasis!')
        except:
            print('Tomographic properties do not hold!')
        finally:
            print('Tomographic properties hold!')
    
    def Covariance(self):
        try:
            for k in F:
                for l in F:
                    for a in F:
                        for b in F:
                            disp = D(k, l)
                            disp_op = disp @ Wootters(a, b) @ disp.conj().T
                            if not np.all(np.isclose(disp_op, Wootters(a+k, b+l))):
                                raise Exception('Not covariant', a, b, k, l)
        finally:
            print('Kernel is covariant!')

    # Testing tomographic probabilities

    def TestProb(self, state):
        w = WW(state, approx=False)
        for mu in F:
            for nu in F:
                wigner_prob = Prob(w, lambda t: mu * t + nu)
                trans_prob  = abs(
                    state.conj() @ mubs306[toInt(mu)+1][:,toInt(nu)]
                )**2
                if not np.isclose(wigner_prob, trans_prob):
                    raise Exception(
                        'Probabilities do not match!',
                        mu, nu, wigner_prob, trans_prob
                    )
        print('Probabilities are correct!')

# https://bunkrr.su/a/wQOwxgpn