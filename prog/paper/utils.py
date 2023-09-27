from sage.all import *
import numpy as np

d = 2**3 # hardcoded sorry

# def checkON(A):
#     return A.conjugate_transpose() * A == identity_matrix(d)

# def checkMUB(A, B):
#     m = (A.conjugate_transpose() * B).apply_map(lambda t: abs(t))
#     return m == (1/sqrt(d) * ones_matrix(d))

# def checkMUBs(mubs):
#     for A in mubs:
#         for B in mubs:
#             if A == B:
#                 if not checkON(A):
#                     raise Exception('Found a basis that is not ON!')
#             else:
#                 if not checkMUB(A, B):
#                     raise Exception('Found two basis that are not MUB!')
#     print('All is good!')

def checkMUBs(mubs):
    d = mubs.shape[1]
    n = int(mubs.shape[0] / d) - 1
    for i in range(n+1):
        for j in range(n+1):
            for k in range(d):
                for l in range(d):
                    v1 = mubs[i*d:(i+1)*d,k]
                    v2 = mubs[j*d:(j+1)*d,l]
                    inp = np.abs(np.conjugate(v1) @ v2)
                    if i == j:
                        if k == l:
                            if not np.isclose(inp, 1):
                                print(inp, i, j, k, l)
                                raise Exception('Not normalized.')
                        else:
                            if not np.isclose(inp, 0):
                                print(inp, i, j, k, l)
                                raise Exception('Not orthogonal.')
                    else:
                        if not np.isclose(inp, 1/np.sqrt(d)):
                            print(inp, i, j, k, l)
                            raise Exception('Not mutually unbiased.')
        print('MUB {} check.'.format(i))
    print('MUBs are good!')

def saveMUBs(mubs, name='mubs.npy'):
    # M = mubs[0].numpy(dtype='complex128')
    M = mubs[0]
    for m in mubs[1:]:
        # M = np.concatenate((M, m.numpy(dtype='complex128')))
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
            t = (O1 @ O2).trace()
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

# For SageMath

# Orthonormality
def isOrthonormal(m):
    return m.conjugate_transpose() * m == Id

# MUBs
def isMUB(m1, m2):
    m = (
        m1.conjugate_transpose() * m2
    ).apply_map(lambda t: abs(t)^2)
    return m == ones_matrix(2^N, 2^N) / (2^N)

def checkMUBs(mubs):
    if type(mubs) != list:
        mubs_list = []
        for i in range(2^N+1):
            mubs_list.append(mubs[(i*2^N):(i+1)*2^N,:])
    else:
        mubs_list = mubs
            
    for i in range(2^N+1):
        for j in range(2^N+1):
            if i == j:
                if not isOrthonormal(mubs_list[i]):
                    raise Exception(
                        'Encountered a basis that is not orthonormal!',
                        i
                    )
            else:
                if not isMUB(mubs_list[i], mubs_list[j]):
                    raise Exception(
                        'Encountered non-MUB pairs of basis!',
                        i, j
                    )
    return True