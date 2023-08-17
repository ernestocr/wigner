from sage.all import *
import numpy as np

p = 2
n = 2
q = p**n

# GR(4,2)
R  = PolynomialRing(Integers(4), 't')
t  = R.gen()
GR = R.quotient(t**2+t+1, 'w')
w  = GR.gen()

def ringTrace(a, b):
    s = GR(0)
    for k in range(n):
        s += a**(2**k) + GR(2)*b**(2**k)
    return s 

# Teichmüller
T = [GR(0)] + [w**k for k in range((2**n - 2) + 1)]

# Really brute force
def brute2adic(x):
    for a in T:
        for b in T:
            if x == a + GR(2)*b:
                return (a,b)

def createMubs(q):
    mubs = [identity_matrix(q)]
    for m in T:
        Bm = zero_matrix(SR, q, q)
        for j, v in enumerate(T):
            for i, w in enumerate(T):
                a, b = brute2adic(w**2 * m + 2 * w * v)
                pwr  = Integer(lift(ringTrace(a, b)))
                Bm[i,j] = 1/sqrt(q) * I**pwr
        mubs.append(Bm)
    return mubs

