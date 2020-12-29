# define atom and module

struct Atom
    name::String
    Z::Int
    Ne::Int
    orbitals::Vector{STONG}
    R::Vec3{Float64}
    Atom(name::AbstractString, Z::Integer, Ne::Integer, orbitals::Vector{STONG{N}}, R::Vec3) where N = new(
        name, Z%Int, Ne%Int, orbitals, float(R)
    )
end

Atom(name::AbstractString, R::Vec3) = begin
    name = titlecase(name)
    number = atomic_number[name]
    Rf = float(R)
    number <= 2 && return Atom(name, number, number, [
        STONG(coefficient_data[name]["1s"], Rf)
    ], Rf)
    number <= 10 && return Atom(name, number, number, [
        STONG(coefficient_data[name]["1s"], Rf),
        STONG(coefficient_data[name]["2s"], Rf),
        STONG(coefficient_data[name]["2p"], Rf; type=:px),
        STONG(coefficient_data[name]["2p"], Rf; type=:py),
        STONG(coefficient_data[name]["2p"], Rf; type=:pz)
    ], Rf)
    error("atomic_number > 10 atoms are not support yet")
end

const atomic_number = Dict{String, Int}(
    "H"  => 1,
    "He" => 2,
    "Li" => 3,
    "Be" => 4,
    "B"  => 5,
    "C"  => 6,
    "N"  => 7,
    "O"  => 8,
    "F"  => 9,
    "Ne" => 10,
    "Na" => 11,
    "Mg" => 12,
    "Al" => 13,
    "Si" => 14,
    "P"  => 15,
    "S"  => 16,
    "Cl" => 17,
    "Ar" => 18,
    "K"  => 19,
    "Ca" => 20
)

hydrogen_atom(R::Vec3) = Atom("H", R)
helium_atom(R::Vec3)   = Atom("He", R)
carbon_atom(R::Vec3)   = Atom("C", R)
oxygen_atom(R::Vec3)   = Atom("O", R)
fluorine_atom(R::Vec3) = Atom("F", R)

struct Molecule
    name::String
    charge::Int
    atoms::Vector{Atom}
    Molecule(name::AbstractString, charge::Integer, atoms::Vector{Atom}) = new(name, charge%Int, atoms)
end

Molecule(name::AbstractString, atoms::Vector{Atom}) = Molecule(name, 0, atoms)
Molecule(atoms::Vector{Atom}) = Molecule("", 0, atoms)