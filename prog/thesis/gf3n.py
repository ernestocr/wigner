from sage.all import *

p = 5 # odd prime
n = 3 # odd
s = 1 # relatively prime to n, 1 <= s < n/2

F = GF(p**n, 'x')
x = F.gen()

FF = [F(0)] + [x**j for j in range(0, p**n - 1)]

l_sols = []
r_sols = []

for k in FF:
    for j in FF:
        if k**2 - j**2 == 0:
            l_sols.append((k, j))
        if k**(p**s+1) - j**(p**s+1) == 0:
            r_sols.append((k, j))

print(l_sols == r_sols)

#sols = 0
#for k in FF:
#    if k**2 == 1:
#        print(k)
#        sols += 1
#print(sols)
