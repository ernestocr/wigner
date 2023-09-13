from sage.all import *
import numpy as np
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import axes3d
import matplotlib.colors as colors
import matplotlib.cm as cm

from utils import TestProbs

class Wigner():

    state  = None
    mubs   = None
    cache  = None
    curves = None

    def __init__(self, d, mubs=None, state=None):
        self.d = d
        self.F = GF(d, 'x')

        if mubs is not None:
            self.LoadMubs(mubs)
        if state is not None:
            self.LoadState(state)

    def __repr__(self):
        d = self.d
        s = (self.state is not None)
        c = (self.cache is not None)
        return f"<WignerFunc dim:{d}, state:{s}, cached:{c}>"

    def CheckState(self, rho):
        if not np.isclose(rho.trace(), 1, rtol=1e-07):
            raise Exception('Density matrix is not trace class.')
        if not np.array_equal(np.conj(rho).T, rho):
            raise Exception('Density matrix is not self-adjoint.')
        return True
    
    def CheckMubs(self, mubs):
        return np.shape(mubs) == ((self.d+1)*self.d, self.d)

    def LoadState(self, rho):
        if not self.CheckState(rho):
            raise Exception('Matrix not a valid state!')
        self.state = rho

        if self.mubs is not None:
            self.WignerMatrix()

    def LoadMubs(self, mubs):
        if not self.CheckMubs(mubs):
            raise Exception('MUBs are not valid!')
        self.mubs = mubs
    
    def LoadCurves(self, curves):
        self.curves = curves

    def Idx(self, k):
        return list(self.F).index(k)
    
    def Delta(self, a, b):
        return float(a == b)
    
    def Proj(self, v):
        return np.outer(v, np.conjugate(v))

    def PhasePointOperator(self, a, b):
        # Hard coded for 306 (standard mubs)
        # Fourier basis first (X op)
        op = self.Proj(self.mubs[:self.d][:, self.Idx(a)])
        for i, k in enumerate(self.F):
            for j, l in enumerate(self.F):
                v = self.mubs[(i+1)*self.d:(i+2)*self.d][:,j]
                op += (self.Delta(b, a * k + l) * self.Proj(v))
        op -= np.eye(self.d)
        return op
    
    def Kernel(self, a, b):
        op = self.Proj(self.mubs[:self.d][:, self.Idx(a)])
        for l, curve in enumerate(self.curves):
            for j, k in enumerate(self.F):
                for t in self.F:
                    d1 = self.Delta(a, t)
                    d2 = self.Delta(b, curve(t) + k)
                    v = self.mubs[(l+1)*self.d:(l+2)*self.d][:,j]
                    op += d1 * d2 * self.Proj(v)
        op -= np.eye(self.d)
        return op

    def Wigner(self, a, b):
        return (
            # self.state @ self.PhasePointOperator(a, b)
            self.state @ self.Kernel(a, b)
        ).trace() / self.d

    def WignerMatrix(self, recalc=False):
        if self.cache is None or recalc:
            m = np.zeros((self.d, self.d))
            for i, a in enumerate(self.F):
                for j, b in enumerate(self.F):
                    m[i,j] = np.real(self.Wigner(a, b))
            # Cache matrix in object
            # self.cache = np.rot90(np.rot90(m, 1).T)
            self.cache = m
        return self.cache
    
    def WignerProb(self, curve):
        s = 0
        w = self.WignerMatrix()

        # ray i <-> mub i+1, curve ij <-> mub i+1, column j ?
        for i, a in enumerate(F):
            for j, b in enumerate(F):
                # delta brute force
                if b == curve(a):
                    s += w[i,j]
        return s


def PlotWignerFunction(w):
    data_array = np.rot90(w.cache)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    x_data, y_data = np.meshgrid(
        np.arange(data_array.shape[1]),
        np.arange(data_array.shape[0])
    )

    x_data = x_data.flatten()
    y_data = y_data.flatten()
    z_data = data_array.flatten()

    ax.set_zlim(min(0, z_data.min()) - 0.05, z_data.max() + 0.05)

    # offset = z_data + np.abs(z_data.min())
    # fracs  = offset.astype(float) / offset.max()
    # norm   = colors.Normalize(fracs.min(), fracs.max())
    # color_vals = cm.viridis(norm(fracs.tolist()))

    ax.bar3d(
        x_data,
        y_data,
        np.zeros(len(z_data)),
        0.9, 0.9, z_data,
        # color = color_vals
    )

    return (fig, ax)

# - - - - -- - - - - - - - - - - - - - -- - - - -- - -
# TESTING

# F = GF(8, 'x')
# x = F.gen()

# # 306
# mubs = np.load('mubs-306.npy')
# w = Wigner(8, mubs)
# w.LoadCurves([
#     lambda t: 0,
#     lambda t: x * t,
#     lambda t: x**2 * t,
#     lambda t: x**3 * t,
#     lambda t: x**4 * t,
#     lambda t: x**5 * t,
#     lambda t: x**6 * t,
#     lambda t: t
# ])

# 162
# mubs = np.load('mubs-162.npy')
# w = Wigner(8, mubs)
# w.LoadCurves([
#     lambda t: x * t**2 + x * t**4,
#     lambda t: x * t + x * t**2 + x * t**4,
#     lambda t: x**2 * t + x * t**2 + x * t**4,
#     lambda t: x**3 * t + x * t**2 + x * t**4,
#     lambda t: x**4 * t + x * t**2 + x * t**4,
#     lambda t: x**5 * t + x * t**2 + x * t**4,
#     lambda t: x**6 * t + x * t**2 + x * t**4,
#     lambda t: t + x * t**2 + x * t**4
# ])

# 234
# mubs = np.load('mubs-234.npy')
# w = Wigner(8, mubs)
# w.LoadCurves([
#     lambda t: 0,
#     lambda t: x**6 * t + x**3 * t**2 + x**5 * t**4,
#     lambda t: x**2 * t + x**5 * t**2 + x**6 * t**4,
#     lambda t: x**4 * t + x**3 * t**2 + x**5 * t**4,
#     lambda t: x**3 * t,
#     lambda t: x**5 * t + x**5 * t**2 + x**6 * t**4,
#     lambda t: x**1 * t + x**2 * t**2 + x**1 * t**4,
#     lambda t: t + x**2 * t**2 + x * t**4
# ])

# from utils import checkPhasePointOperators
# ops = []
# for a in w.F:
#     for b in w.F:
#         # ops.append(w.PhasePointOperator(a,b))
#         ops.append(w.Kernel(a,b))
# checkPhasePointOperators(ops)

# Stabilizer states
# for i in range(9):
#     # w.LoadState(w.Proj(w.mubs[8*i:8*(i+1),0]))
#     w.LoadState(w.Proj(mubs[8*i:8*(i+1),0]))
#     w.WignerMatrix(recalc=True)
#     fig, ax = PlotWignerFunction(w) 
#     plt.show()

# GHZ
# s = np.array([1/np.sqrt(2),0,0,0,0,0,0,1/np.sqrt(2)])
# w.LoadState(w.Proj(s))
# w.WignerMatrix(recalc=True)
# fig, ax = PlotWignerFunction(w) 
# plt.show()

# Verifying the transition probabilities

# s = mubs[32:40][:,5]
# s = np.array([1/np.sqrt(2),0,0,0,0,0,0,1/np.sqrt(2)])
# w.LoadState(w.Proj(s))
# w.WignerMatrix(recalc=True)

# fig, ax = PlotWignerFunction(w)
# plt.show()

# TestProbs(w, w.Proj(mubs[8*1:8*(1+1),0]))