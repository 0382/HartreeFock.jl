# HartreeFock.jl

A Hartree-Fock calculation demo for quantum chemistry.

## Install

You should install `boost, eigen3, libint2` first. Take Ubuntu for example
```bash
sudo apt install libbboost-dev libeigen3-dev libint2-dev
```

Then install this package with `Pkg`
```julia
using Pkg
Pkg.add("https://github.com/0382/HartreeFock.jl.git")
Pkg.build("HartreeFock")
```

## Example

Here is a example to calculate H2O molecule energy
```julia
function H2O(r::Real, θ::Real)
    Molecule("H2O",
        [
            Atom("O", vec3(0, 0, 0)),
            Atom("H", vec3(r*cosd(θ/2),  r*sind(θ/2), 0)),
            Atom("H", vec3(r*cosd(θ/2), -r*sind(θ/2), 0))
        ]
    )
end

println(hartree_fock(H2O(0.99, 105)))
```

For more examples, please see test folder.