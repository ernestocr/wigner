# Testing Gauss sums

from sage.all import *

p = 3
n = 3
d = p**n

F = GF(d, 'x')
x = F.gen()
FF = [F(0)] + [x**j for j in range(27-1)]

omega = exp(2 * pi * I / p)

def char1(m, w, u):
    pwr = (1/F(2)) * m * (w**2 - u**2)
    return omega**int(pwr.trace())

def char2(m, w, u):
    pwr = (1/F(2)) * m * (w**10 - u**10 + m**2 * (w**4 - u**4))
    return omega**int(pwr.trace())

def charSum(w, u, Char):
    return sum([Char(m, w, u) for m in FF])

""" Testing """

#print(charSum(FF[1], FF[14], char1))
#print(charSum(FF[1], FF[14], char2))

#for w in FF:
#    for u in FF:
#        s1 = charSum(w, u, char1)
#        s2 = charSum(w, u, char2)
#        if s1 != s2:
#            raise Exception('Not the same!')
#        print(s1)


