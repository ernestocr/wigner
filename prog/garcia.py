# Studying Garcia and Klimov's 2010 paper

# Sp_{2n}(Z_2) order for qubit systems
def SymplecticGroupSize(n):
    prod = 1
    for k in range(n):
        prod *= (2**(2*(k+1)) - 1) 
    return 2**(n**2) * prod

def MUBsSize(n):
    prod = 1
    for k in range(n):
        prod *= (2**(k+1) + 1)
    return prod
