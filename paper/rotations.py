from sage.all import *
import numpy as np

"""
This code implements the Wigner function
following Muñoz's paper. It is implemented in two
equivalent ways, using the displacement operators 
and using Wootters' kernel.
"""

class Wigner():

    mubs = []

    def __init__(self, N, basis=None, init_mubs=True):
        self.N = N
        self.d = 2**N
        self.F = GF(self.d, 'x')
        self.x = self.F.gen()
        
        self.Id = np.eye(self.d)
        self.FF = self.Fourier()

        # hard coded field basis
        N2BASIS = [self.x, self.x + 1]
        N3BASIS = [self.x**3, self.x**5, self.x**6]

        if not basis:
            if N == 2:
                self.basis = N2BASIS
            elif N == 3:
                self.basis = N3BASIS
            else:
                raise Exception(
                    'Must set basis for dimension {}!'.format(N)
                )
        else:
            self.basis = basis

        if init_mubs:
            self.createMUBs()

    def Proj(self, u, v=None):
        if not v:
            v = u
        return np.outer(u, v.conj().T)

    def toInt(self, k):
        return list(self.F).index(k)

    def chi(self, k):
        return np.exp(np.pi * 1j * int(k.trace()))

    def Fourier(self):
        s = np.zeros((self.d, self.d), dtype='complex128')
        for i, a in enumerate(self.F):
            for j, b in enumerate(self.F):
                s[i,j] = self.chi(a * b) / np.sqrt(self.d)
        return s

    def Z(self, a):
        return np.diag([self.chi(a * k) for k in self.F])
    
    def X(self, b):
        return self.FF.conj().T @ self.Z(b) @ self. FF

    def components(self, k):
        return [(k * m).trace() for m in self.basis]

    # Bit counting function
    def h(self, k):
        return sum([int(bit) for bit in self.components(k)])

    # Rotation coefficients
    def c(self, alpha, xi, p=1):
        return (-1j)**(self.h(alpha**p * np.sqrt(xi)**p))

    # Phase defined by the rotation coeffs
    def phi(self, tau, nu, p=1):
        t = type(self.F(0))
        if type(tau) != t:
            tau = self.F(tau)
        if type(nu) != t:
            nu = self.F(nu)
            
        if tau == self.F(0):
            return 1
        return self.c(tau, tau**-1 * nu, p)

    # Displacement operators with phase
    def D(self, a, b, p=1):
        return self.phi(a, b, p) * self.Z(a) @ self.X(b)

    # Rotation operators
    def V(self, xi, p=1):
        s = np.zeros((self.d, self.d), dtype='complex128')
        for i, k in enumerate(self.F):
            s += self.c(k, xi, p) * self.Proj(self.FF[:,i])
        return s

    # Displacement kernel
    def Delta(self, a, b):
        s = np.zeros((self.d, self.d), dtype='complex128')
        for gamma in self.F:
            for delta in self.F:
                char = self.chi(a * delta + b * gamma)
                s += char * self.D(gamma, delta) / self.d
        return s

    def createMUBs(self):
        self.mubs = [self.FF] + [self.V(k, p=1) for k in self.F]

    def Wootters(self, a, b):
        op = self.Proj(self.FF[:, self.toInt(a)])
        for xi in self.F:
            for nu in self.F:
                d = int(b == xi * a + nu)
                v = self.mubs[self.toInt(xi)+1][:, self.toInt(nu)]
                op += d * self.Proj(v)
        return op - self.Id

    def Wigner(self, state, a, b, kernel):
        return (state @ kernel(a, b)).trace()

    def selectKernel(self, ker):
        if ker == 'wootters':
            return self.Wootters
        return self.Delta
    
    def WignerMatrix(self, state, ker='delta'):
        kernel = self.selectKernel(ker)
        
        W = np.zeros((self.d, self.d), dtype='float64')
        for i, a in enumerate(self.F):
            for j, b in enumerate(self.F):
                W[i, j] = real(
                    self.Wigner(state, a, b, kernel)
                ) / self.d
        return W
