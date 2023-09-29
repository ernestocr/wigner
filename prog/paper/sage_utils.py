from sage.all import *

N = 3
Id = identity_matrix(SR, 2^N)

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