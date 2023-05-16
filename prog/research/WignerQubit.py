import numpy as np
import galois
import matplotlib.pyplot as plt

d = 2
n = 2
N = d**n
GF = galois.GF(N)

def expand(v, dual):
    """ Expand field element v with respect to some given basis. """
    return np.array([(x * dual).field_trace() for x in v])

# 2x2 Pauli matrices in standard basis.
X = np.array([[0, 1], [1, 0]]) # unit horizontal translation
Z = np.array([[1, 0], [0, -1]]) # unit vertical translation

# Wootters choice for 2^2: (w,1)
basis = GF([2,1])
dual_base = GF([1,3])
def translation(v):
    """ Translation operator for a given phase-space translation. """
    coeffs = expand(v, dual_base)
    T = [1]
    for k in coeffs.T:
        M1 = np.linalg.matrix_power(X, k[0])
        M2 = np.linalg.matrix_power(Z, k[1])
        T = np.kron(T, M1 @ M2)
    return T

def projector(v):
    """ Rank-1 projection operator: |v><v|."""
    return np.kron(np.conj(v), v.reshape((-1,1)))

# MUBs for two qubits corresponding to translations,
# that Wootters used.  We don't know what method he
# used to diagonalize, if he diagonalized at all.

# (0,1)
s1 = np.eye(N)

# (1,0)
s2 = np.array(
    [[1,1,1,1],
     [1,-1,1,-1],
     [1,1,-1,-1],
     [1,-1,-1,1]]
).T / 2

# (1,1)
s3 = np.array(
    [[1,-1j, 1j, 1],
     [1,1j,1j,-1],
     [1,-1j,-1j,-1],
     [1,1j,-1j,1]]
).T / 2

# (1,w)
s4 = np.array(
    [[1,1,1j,-1j],
     [1,-1,1j,1j],
     [1,1,-1j,1j],
     [1,-1,-1j,-1j]]
).T / 2

# (1,w+1)
s5 = np.array(
    [[1,-1j,1,1j],
     [1,1j,1,-1j],
     [1,-1j,-1,-1j],
     [1,1j,-1,1j]]
).T / 2

striations = [s1, s2, s3, s4, s5]

def isInvariant(line, t):
    """ Verify if line is invariant under translation t. """
    l1 = set(map(tuple, np.array(line)))
    l2 = set(map(tuple, np.array(line + t)))
    return l1 == l2

def getRay(t):
    """
    The ray consists of elements (sq,sp) for some (q,p),
    and all elements s of GF.
    """
    return GF([s * t for s in GF.elements])

def quantumNet(line):
    """ Returns the rank-1 projector associated to the line. """
    # first the find striation where the line belongs to,
    # i.e., which basic translation leaves the line invariant.
    ts = [[0,1], [1,0], [1,1], [1,2], [1,3]]
    idx = None
    for (i,t) in enumerate(ts):
        if isInvariant(line, GF(t)):
            idx = i
            break
    # get the ray of the striation
    ray = getRay(GF(ts[idx]))
    # calculate the displacement from ray to line
    T = translation(GF(line[0] - ray[0]))
    # get first basis element of the corresponding striation (ray)
    basis_1 = striations[idx][:,0]
    # translate the ray projector according to the displacement T
    return T @ projector(basis_1) @ np.conj(T).T

def OriginPointOperator():
    """
    Calculates A(0,0) which is just the sum of the projectors
    associated to the rays, which in turn are just the
    projectors of the first element of our MUBs.
    """
    # ts = [[0,1], [1,0], [1,1], [1,2], [1,3]]
    Op = np.zeros((N,N))
    # for (i,t) in enumerate(ts):
    for i in range(N+1):
        # Op = Op + quantumNet(getRay(GF(t)))
        Op = Op + projector(striations[i][:,0])
    return (Op - np.eye(N)) / N

def PointOperator(v):
    """
    Returns the phase-point operator corresponding to v
    by translating the operator A(0,0).
    """
    T = translation(GF(v))
    return T @ OriginPointOperator() @ np.conj(T).T

def Wigner(rho):
    """
    Returns the Wigner function (as a matrix) for
    a density matrix rho.
    $$
    W(i,j) = Tr(rho @ A(i,j)).
    $$
    """
    W = np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            W[i,j] = np.real((PointOperator([i,j]) @ rho).trace())
    # rotate matrix to match cartesian coordinates.
    return np.rot90(W)

# -------------------------------------------------------------------
# EXAMPLES.

print("|00>")
v = np.kron([1,0], [1,0])
print(Wigner(projector(v)))

print("|11>")
v = np.kron([0,1], [0,1])
print(Wigner(projector(v)))

print("|1->>")
v = np.kron([1,0], [1,1]/np.sqrt(2))
print(Wigner(projector(v)))

print("1/sqrt 2 |10> - |01>")
v = (np.kron([1,0], [0,1]) - np.kron([0,1],[1,0])) / np.sqrt(2)
print(Wigner(projector(v)))

print("Bell state +")
bell2 = (np.kron([1,0], [1,0]) + np.kron([0,1],[0,1])) / np.sqrt(2)
print(Wigner(projector(bell2)))