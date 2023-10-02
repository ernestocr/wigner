from sage.all import *
import numpy as np

def testMubs(mubs):
    d = mubs[0].nrows()
    for b, B in enumerate(mubs):
        print('Testing set', b)
        for v, V in enumerate(mubs):
            if b != v:
                for j in range(d):
                    for k in range(d):
                        r = (transpose(B[:,j].conjugate()) * V[:,k])[0,0]
                        r = (r.abs()**2).simplify()
                        if r != 1/d:
                            print(r, b, v, j, k)
                            raise Exception('Not MUBs!')
    return True

def testMubs2(mubs):
    d = mubs[0].nrows()
    y = (1/d) * ones_matrix(SR, d)
    for b, B in enumerate(mubs):
        print('Testing set', b)
        for v, V in enumerate(mubs):
            if b != v:
                x = B.conjugate_transpose() * V
                x = x.apply_map(lambda z: (z.abs()**2).simplify())
                if x != y:
                    print(b, v)
                    raise Exception('Not MUBs!')
    return True

# def testMubsNp(mubs):
#     d = mubs[0].nrows()
#     mubs = [m.numpy(dtype='complex64) for m in mubs]
#     cmp  = np.ones((d,d)) / d
#     for b, B in enumerate(mubs):
#         for v, V in enumerate(mubs):
#             if b != v:
#                 X = B.conj().T @ V
#                 X = np.abs(X)**2
#                 if not np.allclose(X, cmp):
#                     print(b, v)
#                     raise Exception('Not MUBs!')
#     return True

def saveMubs(mubs, path):
    nmubs = [m.numpy(dtype='complex64') for m in mubs]
    np.save(path, np.concatenate(nmubs))
    print('Saved to', path)
    pass

def proj(v):
    d = len(v)
    v = v.reshape((d,1))
    return np.kron(v, v.conj().T)
