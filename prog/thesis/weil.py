from sage.all import *

p = 3
n = 3
d = p**n

F = GF(d, 'x')
x = F.gen()
FF = [F(0)] + [x**j for j in range(27-1)]

omega = exp(2 * pi * I / p)

# The Weil sum for f(x) = 1/2 x ( w^2 - u^2 ) is 0
# whenever w^2 =\= u^2, since f(x) is linear.
# So the roots of the equation w^2 - u^2 gives us the
# value of 27.

sols_left  = []
sols_right = []
for i, w in enumerate(FF):
    for j, u in enumerate(FF):
        if w**2 - u**2 == 0:
            sols_left.append((i, j))
        #if (w**4 - u**4)**3 + F(2)**9*(w**10 - u**10)**9 == 0:
        if w**4 - u**4 == 0:
            sols_right.append((i, j))

print('Number of solutions for left side: ', len(sols_left))
print('Number of solutions for right side: ', len(sols_right))
print(sols_left == sols_right)
