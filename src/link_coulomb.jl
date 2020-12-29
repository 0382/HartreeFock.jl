module Coulomb
    using CxxWrap
    @wrapmodule(joinpath(@__DIR__, "../deps/libcoulomb.so"))
    function __init__()
        @initcxx
        initialize()
    end
end

import .Coulomb

Gaussian_size(::Gaussian_s) = 1
Gaussian_size(::Gaussian_p) = 3
Gaussian_size(s::STO3G) = Gaussian_size(s.g[1])

function make_shells(basis::Vector{STO3G})
    l = Vector{Cint}()
    alpha = Vector()
    d = Vector()
    R = Vector()
    i = 1
    Nshell = 0
    while i <= length(basis)
        sto = basis[i]
        g = sto.g[1]
        step = Gaussian_size(sto)
        push!(l, step == 1 ? 0 : 1)
        push!(alpha, [sto.g[1].α, sto.g[2].α, sto.g[3].α])
        push!(d, [sto.d...])
        push!(R, [g.R[1], g.R[2], g.R[3]])
        i += step
        Nshell += 1
    end
    return Nshell, l, vcat(alpha...), vcat(d...), vcat(R...)
end

function Ve!(ve::Array, basis::Vector{STO3G})
    Nshell, l, alpha, d, R = make_shells(basis)
    Coulomb.Ve(
        pointer(ve), Nshell, pointer(l), pointer(alpha), pointer(d), pointer(R)
    )
end

function Veff!(veff::Array, ρ::Array, Nshell::Int, l::Vector{Cint}, alpha::Vector{Float64}, d::Vector{Float64}, R::Vector{Float64})
    Coulomb.Veff(
        pointer(veff), pointer(ρ), Nshell, pointer(l), pointer(alpha), pointer(d), pointer(R)
    )
end