from sage.all import *
import numpy as np
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import axes3d
import matplotlib.colors as colors
import matplotlib.cm as cm

class Wigner():

    state = None
    mubs  = None
    cache = None

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
        # Add all operators whose lines contain (a,b)
        for i, k in enumerate(self.F):
            for j, l in enumerate(self.F):
                v = self.mubs[(i+1)*self.d:(i+2)*self.d][:,j]
                op += (self.Delta(b, a * k + l) * self.Proj(v))
        op -= np.eye(self.d)
        return op
    
    def Kernel(self, a, b):
        pass

    def Wigner(self, a, b):
        return (
            self.state @ self.PhasePointOperator(a, b)
        ).trace() / self.d

    def WignerMatrix(self, recalc=False):
        if self.cache is None or recalc:
            m = np.zeros((self.d, self.d))
            for i, a in enumerate(self.F):
                for j, b in enumerate(self.F):
                    m[i,j] = np.real(self.Wigner(a, b))
            # Cache matrix in object
            self.cache = np.rot90(np.rot90(m, 1).T)
        return self.cache

def PlotWignerFunction(w):
    data_array = w.cache

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

mubs = np.load('mubs-306.npy')
w = Wigner(8, mubs)

# from utils import checkPhasePointOperators

# ops = []
# for a in w.F:
#     for b in w.F:
#         ops.append(w.PhasePointOperator(a,b))
# checkPhasePointOperators(ops)

# i = 8
# w.LoadState(w.Proj(w.mubs[8*i:8*(i+1),0]))
s = np.array([1/np.sqrt(2),0,0,0,0,0,0,1/np.sqrt(2)])
w.LoadState(w.Proj(s))
fig, ax = PlotWignerFunction(w) 
plt.show()