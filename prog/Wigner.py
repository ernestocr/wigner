from sage.all import *
import numpy as np

class Wigner():

    def __init__(self, field, mubs, op=None):
        self.field = field
        self.mubs = mubs

        self.op = self.Desargues
        if op == 'kantor':
            self.op = self.Kantor
        if op == 'albert':
            self.op = self.Albert

    def toInt(self, e):
        return list(self.field).index(e)
    
    def npProj(self, v):
        d = len(v)
        return np.kron(v.reshape((d,1)), v.conj())

    # Presemifield operation
    def Desargues(self, m, u):
        return m*u
    
    def Kantor(self, m, u):
        return m**2*u + m*u.trace() + (m*u).trace()
    
    def Albert(self, m, u):
        return m*u**9 + m**3*u**3
    
    # Geometry
    def Spread(self):
        lines = [[(self.field(0), u) for u in self.field]]
        for m in self.field:
            lines.append([(u, self.op(m, u)) for u in self.field])
        return lines

    def PointLines(self, x, y):
        bs = [y - self.op(m, x) for m in self.field]
        
        lines = [[(x, u) for u in self.field]]
        for i, m in enumerate(self.field):
            line = [(x, self.op(m, x) + bs[i]) for x in self.field]
            lines.append(line)

        return lines
    
    # Wigner function
    def A(self, x, y):
        d = self.field.order()
        bs = [x] + [y - self.op(m, x) for m in self.field]
        bi = [self.toInt(b) for b in bs]
        
        op = np.zeros((d,d), dtype='complex64')
        for k in range(d + 1):
            B = self.mubs[d*k:d*(k+1), :]
            p = self.npProj(B[:, bi[k]])
            op += p
            
        return op - np.eye(d)
    
    def Wigner(self, rho, x, y):
        d = self.field.order()
        return (rho @ self.A(x, y)).trace() / d

    def WignerMatrix(self, rho):
        d = rho.shape[0]
        W = np.zeros((d,d))
        for i, x in enumerate(self.field):
            for j, y in enumerate(self.field):
                W[i,j] = np.real(self.Wigner(rho, x, y))
        return np.rot90(W)