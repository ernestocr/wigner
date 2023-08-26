from sage.all import *

d = 2**3 # hardcoded sorry

def checkON(A):
    return A.conjugate_transpose() * A == identity_matrix(d)

def checkMUB(A, B):
    m = (A.conjugate_transpose() * B).apply_map(lambda t: abs(t))
    return m == (1/sqrt(d) * ones_matrix(d))

def checkMUBs(mubs):
    for A in mubs:
        for B in mubs:
            if A == B:
                if not checkON(A):
                    raise Exception('Found a basis that is not ON!')
            else:
                if not checkMUB(A, B):
                    raise Exception('Found two basis that are not MUB!')
    print('All is good!')