include("../src/HartreeFock.jl")

using .HartreeFock
using BenchmarkTools

# define some unit vector
const n = (
    vec3(  1.,    0., 0.),
    vec3( 1/2,  √3/2, 0.),
    vec3(-1/2,  √3/2, 0.),
    vec3( -1.,    0., 0.),
    vec3(-1/2, -√3/2, 0.),
    vec3( 1/2, -√3/2, 0.)
)

# dCC: bond length between adjacent `C` Atom
# dCH: bond length between `C` and `H`
function Benzene(dCC::Float64, dCH::Float64)::Molecule
    Molecule("C6H6", vcat(
        collect(ntuple(i->Atom("C", dCC * n[i]), 6)),
        collect(ntuple(i->Atom("H", (dCC+dCH) * n[i]), 6))
    ))
end

function Benzene_pes(dCC_range, dCH_range)
    NCC = length(dCC_range)
    NCH = length(dCH_range)
    Emat = zeros(Float64, NCC, NCH)
    Threads.@threads for i in 1:NCC
    # for i in 1:NCC
        for j in 1:NCH
            benzene = Benzene(dCC_range[i] / a0, dCH_range[j] / a0)
            total_energy = hartree_fock(benzene)
            Emat[i, j] = total_energy
        end
    end
    i, j = argmin(Emat).I
    println("minimal energy is $(Emat[i,j]), with dCC = $(dCC_range[i])Å, dCH = $(dCH_range[j])Å")
end

#display(@benchmark Benzene_pes(1.3:0.02:1.5, 1:0.02:1.2))
Benzene_pes(1.3:0.02:1.5, 1:0.02:1.2)