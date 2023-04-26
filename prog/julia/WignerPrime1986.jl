"""
    Wigner(ρ)

Compute the discrete Wigner function for a qudit in
the state ρ, for prime state space dimension `d`,
following Wootters' original method (1986).

# Arguments
- `ρ::Matrix`: the density operator of the state.

# Examples.

(+1)-eigenstate of σz
```julia-repl
julia> v  = [1; 0]
julia> vv = kron(adjoint(v), v)
julia> Wigner(vv)
```

(+1)-eigenstate of (1/√2)(σz + σx)
```julia-repl
julia> v  = eigvecs((σz+σx)/√2)[:,2]
julia> vv = kron(adjoint(v), v)
julia> Wigner(vv)
```
"""

using AbstractAlgebra
using LinearAlgebra
using Primes
using Plots
using LaTeXStrings

const TOL = 1e-7 # for verification of density matrix.

# Pauli matrices in the standard basis (for N = 2)
σx = [0 1; 1 0]
σy = [0 -im; im 0]
σz = [1 0; 0 -1]

# Kronecker delta 
δ(x, y) = ==(x,y)

function PointOperator(
    a₁::AbstractAlgebra.GFElem{Int64},
    a₂::AbstractAlgebra.GFElem{Int64},
    F::AbstractAlgebra.GFField{Int64})
    
    N = length(F)

    if N == 2
        return (
            (-1)^lift(a₁) * σz +
            (-1)^lift(a₂) * σx +
            (-1)^lift(a₁ + a₂) * σy +
            I
        ) / 2 
    end
    
    A = zeros(Complex{Float64}, N, N)
    for (i, k) in enumerate(F)
        for (j, l) in enumerate(F)
            A[i, j] = δ(2*a₁, k+l) * exp((2*π*im / N) * lift(a₂*(k-l)))
        end
    end
    return A
end

function isPositiveSemidefinite(A)
    eigs = eigvals(A)
    return ishermitian(A) && all(>=(-TOL), eigs)
end

function isDensity(A)
    return abs(tr(A) - 1.0) < TOL && isPositiveSemidefinite(A)
end

function Wigner(ρ)
    N = size(ρ)[1]
    
    !isprime(N) && error("N is not a prime number!")
    !isDensity(ρ) && error("Matrix is not a valid density operator!")

    F = GF(N) # Finite field Z_n

    W = zeros(Float64, N, N)
    for (i, k) in enumerate(F)
       for (j, l) in enumerate(F)
            W[i,j] = real(tr(ρ * PointOperator(k, l, F)))
       end
    end

    return convert(Matrix{Float64}, CartesianOrder(W / N))
end

function CartesianOrder(W)
    return W
    # return view(transpose(W), [2,1], :)
end

function heat(W)
    N = size(W)[1]
    heatmap(
        CartesianOrder(W), # reverse order
        xticks = (1:N, 0:N-1),
        yticks = (1:N, 0:N-1),
        xlabel = L"q",
        ylabel = L"p",
        # aspect_ratio = :equal,
        c = :viridis,
        legend = false,
        colorbar = true
    )
    vline!(0.5:(N+0.5), c = :black)
    hline!(0.5:(N+0.5), c = :black)
end

v  = eigvecs((σz+σx)/√2)[:,2]
vv = kron(adjoint(v), v)
heat(Wigner(vv))

heat(Wigner(I(5)/5))

v1 = [0; 0; 0; 1; 0]
v2 = [0; 1; 0; 0; 0]
v  = (v1 + v2) / sqrt(2)
vv = kron(adjoint(v), v)
heat(Wigner(vv))