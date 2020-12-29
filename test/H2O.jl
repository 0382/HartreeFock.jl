include("../src/HartreeFock.jl")

using .HartreeFock
using BenchmarkTools

function H2O(r::Real, θ::Real)
    Molecule("H2O",
        [
            Atom("O", vec3(0, 0, 0)),
            Atom("H", vec3(r*cosd(θ/2),  r*sind(θ/2), 0)),
            Atom("H", vec3(r*cosd(θ/2), -r*sind(θ/2), 0))
        ]
    )
end

function H2O_pes(r_range::AbstractArray, θ_range::AbstractArray)
    Nr = length(r_range)
    Nt = length(θ_range)
    total_energy = zeros(Float64, Nr, Nt)
    Threads.@threads for i in 1:Nr
    # for i in 1:Nr
        r = r_range[i] / a0
        for j in  1:Nt
            θ = θ_range[j]
            e = hartree_fock(H2O(r, θ))
            total_energy[i, j] = e
        end
    end
    i, j = argmin(total_energy).I
    println("minimal energy of H2O is $(total_energy[i,j]*hartree)eV, with the bond length = $(r_range[i])Å, bond angle = $(θ_range[j])ᵒ")
end

# display(@benchmark H2O_pes(0.9:0.01:1.1, 85:115))
H2O_pes(0.9:0.01:1.1, 85:115)