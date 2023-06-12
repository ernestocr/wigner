from sage.all import *
import numpy as np
import matplotlib.pyplot as plt

class Wigner():

    def __init__(self, field, mubs, op=None):
        self.field = field
        self.mubs = mubs

        self.op = self.Desargues
        if op != None:
            if op.lower() == 'kantor':
                self.op = self.Kantor
            if op.lower() == 'albert':
                self.op = self.Albert

    def toInt(self, e):
        return list(self.field).index(e)
    
    def npProj(self, v):
        d = len(v)
        v = v.reshape((d,1))
        return np.kron(v, v.conj().T)

    # Presemifield operation
    def Desargues(self, m, u):
        return m*u
    
    def Kantor(self, m, u):
        return m**2*u + m*u.trace() + (m*u).trace()
    
    def Albert(self, m, u):
        return m*u**9 + m**3*u**3
    
    # Geometry
    def Spread(self):
        lines = [[(self.field[0], u) for u in self.field]]
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
        d = len(self.field)
        bs = [x] + [y - self.op(m, x) for m in self.field]
        bi = [self.toInt(b) for b in bs]
        
        a = np.zeros((d,d), dtype='complex64')
        for k in range(d + 1):
            B = self.mubs[d*k:d*(k+1)]
            p = self.npProj(B[:, bi[k]])
            a += p
            
        return a - np.eye(d)
    
    def Wigner(self, rho, x, y):
        d = len(self.field)
        return (rho @ self.A(x, y)).trace() / d

    def WignerMatrix(self, rho):
        d = rho.shape[0]
        W = np.zeros((d,d))
        for i, x in enumerate(self.field):
            for j, y in enumerate(self.field):
                W[i,j] = np.real(self.Wigner(rho, x, y))
        return np.rot90(W)

def plotHeat(W, ax=None):
    if ax == None:
        fig, ax = plt.subplots()
    
    ax.imshow(np.rot90(W.T), origin='lower')
    return ax

def npProj(v):
    d = len(v)
    v = v.reshape((d,1))
    return np.kron(v, v.conj().T)