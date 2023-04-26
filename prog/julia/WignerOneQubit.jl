"""
Wigner function for one qubit using Gibbons and 
Wootters newer version.
"""

# Bases[BASE][:, BASIS_ELEMENT]
# ==
# Bases[STRIATION][:, LINE]

# An indexing scheme would be much faster. For example if we
# fix the way we list the lines of a striation by slope  and
# intercept in some predetermined order, then we can order 
# the basis elements according to the translations from the ray.
# But this would require calculating the translation of the ray 
# and basis element any way! So it would have to be done one time
# upon instantiation of some Wigner object. And then simply identify
# the projector by the line given its slope and intercept, and 
# proceed to look up on the fixed MUBs according to the indexing.

# A slower but easier and more direct way to obtain the projectors
# would be to simply calculate the translation necessary to 
# obtain the selected line from each ray. And then use the
# translation operators to get the basis element.

using LinearAlgebra
using AbstractAlgebra
using LaTeXStrings
using Plots

F  = GF(2)

Bases = [
    [1 0;
    0 1], # m = Inf
    [1/sqrt(2) 1/sqrt(2); 
    1/sqrt(2) -1/sqrt(2)], # m = 0
    [1/sqrt(2) 1/sqrt(2); 
    im/sqrt(2) -im/sqrt(2)] # m = 1
]

X = [0 1; 1 0]
Z = [1 0; 0 -1]

function Translation(v)
    return X^v[1] * Z^v[2]
end

function Projector(v)
    return kron(adjoint(v), v)
end

# Ray <-> First basis element
function QuantumNet(
    m::AbstractAlgebra.GFElem{Int64},
    c::AbstractAlgebra.GFElem{Int64})
    idx = findfirst(==(m), [e for e in F])
    ray = Projector(Bases[idx + 1][:,1])
    T   = Translation([0,lift(c)])
    return T * ray * adjoint(T)
end

function QuantumNet(m::Float64, c::AbstractAlgebra.GFElem{Int64})
    ray = Projector(Bases[1][:,1])
    T   = Translation([lift(c), 0])
    return T * ray * adjoint(T)
end

function PointOperator(
    q::AbstractAlgebra.GFElem{Int64},
    p::AbstractAlgebra.GFElem{Int64})
    Op = QuantumNet(Inf, q)
    for m in F
        for c in F
            if p == m * q + c
                Op += QuantumNet(m, c)
            end
        end
    end
    return (Op - I(2)) / 2
end

function Wigner(ρ)
    W = zeros(2, 2)
    for (i, k) in enumerate(F)
        for (j, l) in enumerate(F)
            W[i, j] = real(tr(ρ * PointOperator(k, l)))
        end
    end
    return rotl90(W, 1)
end

heatmap(Wigner(I(2)/2))
heatmap(Wigner(Projector([1; 0])))
heatmap(Wigner(Projector(Bases[2][:,2])))
heatmap(Wigner(Projector(Bases[3][:,1])))