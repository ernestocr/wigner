# Following Vourdas, for arbitrary ring Z(d).
# Caveat: formalism only works for d > 2.
# Question: What do we lose since we have not
# restricted to power of prime?

using AbstractAlgebra
using LinearAlgebra

const d = 3
const F = GF(3)

# Using the Fourier transform

function ω(α::AbstractAlgebra.GFElem{Int64})
    return exp(im * (2 * π * lift(α)) / d)
end

function θX(
    m::AbstractAlgebra.GFElem{Int64},
    n::AbstractAlgebra.GFElem{Int64},
    θ)
    return adjoint(I(d)[:,lift(m)+1]) * θ * I(d)[:,lift(n)+1]
end

function W(
    α::AbstractAlgebra.GFElem{Int64},
    β::AbstractAlgebra.GFElem{Int64},
    θ)
    return sum([ω(2*α*n) * θX(β-n, β+n, θ) for n in F])
end

function W(θ)
    w = zeros(Complex{BigFloat}, d, d)
    for (i, α) in enumerate(F)
        for (j, β) in enumerate(F)
            w[i,j] = W(α, β, θ)
        end
    end
    return real(w)
end

# Verbatim

# Fourier transform
FF = zeros(Complex, d, d)
for (i, m) in enumerate(F)
    for (j, n) in enumerate(F)
        FF += ω(m*n) * I(d)[:,i] * adjoint(I(d)[:,j])
    end
end
FF = FF / sqrt(d)

# X and Z
x = sum([(n-1) * I(d)[:,n] * adjoint(I(d)[:,n]) for n in 1:d])
Z = exp(im*2*π / d * x)
X = adjoint(FF) * Z * FF

# Displacement
function D(α::AbstractAlgebra.GFElem, β::AbstractAlgebra.GFElem)
    return ω(-inv(F(2)) * α * β) * X^(lift(α)) * Z^(lift(β))
end

# Parity 
function P(α::AbstractAlgebra.GFElem, β::AbstractAlgebra.GFElem)
    return D(α,β) * FF^2 * adjoint(D(α,β))
end

# W
function W2(α, β, θ)
    return tr(θ * P(α, β))
end

function W2(θ)
    w = zeros(Complex{BigFloat}, d, d)
    for (i, α) in enumerate(F)
        for (j, β) in enumerate(F)
            w[i,j] = W2(α, β, θ)
        end
    end
    return transpose(real(w))
end