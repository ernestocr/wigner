from sage.all import *
import numpy as np

# d = 2**5 # hardcoded sorry

def checkON(A):
    d = A.shape[0]
    return np.allclose(A.conj().T @ A, np.eye(d))

def checkMUB(A, B):
    d = A.shape[0]
    m = np.abs(A.conj().T @ B)**2
    return np.allclose(m, np.ones((d, d)) / d)

def checkMUBs(mubs):
    d = mubs[0].shape[1]
    for i, A in enumerate(mubs):
        for j, B in enumerate(mubs):
            if i == j:
                if not checkON(A):
                    raise Exception('Found a basis that is not ON!', i)
            else:
                if not checkMUB(A, B):
                    raise Exception('Found two basis that are not MUB!', i, j)
    print('MUBs are good!')

def saveMUBs(mubs, name='mubs.npy'):
    M = mubs[0]
    for m in mubs[1:]:
        M = np.concatenate((M, m))
    np.save(name, M)

def checkPhasePointOperators(ops):
    d = ops[0].shape[0]
    for O1 in ops:
        if not np.all(O1.conj().T == O1):
            raise Exception('Operator is not self-adjoint.')
        if not np.isclose(O1.trace(), 1):
            raise Exception('Operator is not unit trace.')
    print('Operators are self-adjoint and unit trace!')

    for i, O1 in enumerate(ops):
        for j, O2 in enumerate(ops):
            t = (O1.conj().T @ O2).trace()
            if np.all(O1 == O2):
                if not np.isclose(t, d):
                    print(t, i, j)
                    raise Exception('Not orthonormal.')
            else:
                if not np.isclose(t, 0):
                    print(t, i, j)
                    raise Exception('Not orthonormal.')
    print('Operators are orthonormal!')

# For verifying probabilities

def TransitionProb(s1, s2):
    return np.abs(np.dot(s1, s2.conjugate()))**2

def TestProbs(w, state, F=None):
    if F == None:
        d = w.d
        F = GF(d, 'x')

    for i, l in enumerate(w.curves):
        for j, k in enumerate(F):
            pl = lambda t: l(t) + k
            wp = w.WignerProb(pl)
            tp = TransitionProb(state, w.mubs[(i+1)*8:(i+2)*8][:,j])
            if not np.isclose(tp, wp):
                print(i, j, tp, wp)
                raise Exception('Probabilities do not match!')

    print('Probabilities match!')
