"Restrict Hartree-Fock"
function hartree_fock(mol::Molecule;
        max_iter::Int = 100, δe::Float64 = 1e-10,
        large_system::Bool = false,
        return_density::Bool = false
    )
    
    if large_system == false
        small_system_hartree_fock(mol; max_iter = max_iter, δe = δe, return_density = return_density)
    else
        large_system_hartree_fock(mol; max_iter = max_iter, δe = δe, return_density = return_density)
    end
end

function small_system_hartree_fock(mol::Molecule;
        max_iter = 100, δe = 1e-10,
        return_density::Bool = false
    )
    basis = Vector{STO3G}()
    Ne::Int = -mol.charge
    for atom::Atom in mol.atoms
        Ne += atom.Ne
        for orbital in atom.orbitals
            push!(basis, orbital)
        end
    end
    @assert iseven(Ne) "odd electron number is not suppport yet"
    Nb::Int = length(basis)
    Nn::Int = length(mol.atoms)
    
    # overlap matrix
    S = zeros(Float64, Nb, Nb) + I
    for i = 1:Nb
        for j = (i+1):Nb
            S[i, j] = S[j,i] = overlap_integral(basis[i], basis[j])
        end
    end
    
    # kinetic matrix
    T = zeros(Float64, Nb, Nb)
    for i = 1:Nb
        for j = i:Nb
            T[i,j] = T[j,i] = kinetic_integral(basis[i], basis[j])
        end
    end

    # nuclear elergy
    nuclear_energy = 0.0
    for I = 1:Nn
        for J = (I+1):Nn
            nuclear_energy += mol.atoms[I].Z * mol.atoms[J].Z / abs(mol.atoms[I].R - mol.atoms[J].R)
        end
    end

    # nuclear attraction matrix
    Vn = zeros(Float64, Nb, Nb)
    for atom in mol.atoms
        for i = 1:Nb
            for j = i:Nb
                _v = nuclear_acctaction_integral(atom.Z, atom.R, basis[i], basis[j])
                Vn[i,j] += _v
                if i != j
                    Vn[j, i] += _v
                end
            end
        end
    end

    # one body matrix
    H0 = T + Vn

    # use libint2 for electron repulsion integral instead
    Ve = zeros(Float64, Nb, Nb, Nb, Nb)
    Ve!(Ve, basis)

    d, U = eigen(Symmetric(S))
    X = U * diagm(d.^(-1/2)) * U'

    ρ = zeros(Float64, Nb, Nb)  # density matrix
    G = zeros(Float64, Nb, Nb)  # effective interaction matrix    
    F = zeros(Float64, Nb, Nb)  # Fock matrix
    H = zeros(Float64, Nb, Nb)  # orthogonal Fock matrix
    C = zeros(Float64, Nb, Nb)  # orthogonal coefficients
    Cx = zeros(Float64, Nb, Nb) # real coefficients

    electron_energy = 0.0
    total_energy = 0.0
    old_energy = Inf64

    # start hartree-fock iteration
    iter_times = 0
    for hf_iter = 1:max_iter
        @. F = H0 + G
        H .= X * F * X
        C .= eigvecs(H)
        Cx .= real(X * C)

        # update density matrix
        ρ .= 0
        for i = 1:div(Ne, 2)
            for α = 1:Nb
                for a = 1:Nb
                    @inbounds ρ[α, a] += 2*Cx[α, i] * Cx[a, i]
                end
            end
        end
        
        # update energy
        electron_energy = 0.0
        for α = 1:Nb
            for a = 1:Nb
                @inbounds electron_energy += ρ[a, α] * (H0[α, a] + F[α, a])
            end
        end
        electron_energy *= 0.5

        if abs(old_energy - electron_energy) < δe
            break
        end

        old_energy = electron_energy

        # update effective interaction matrix
        G .= 0
        for α = 1:Nb
            for β = 1:Nb
                for b = 1:Nb
                    for a = 1:Nb
                        G[α, a] += (Ve[α, β, b, a] - 0.5Ve[α, β, a, b]) * ρ[b, β]
                    end
                end
            end
        end
        iter_times = hf_iter
    end
    
    total_energy = electron_energy + nuclear_energy
    if iter_times == max_iter
        println("Warning: hartree_fock iterator time reaches `max_iter`, but may not reach the convergence condition")
    end
    if return_density
        return total_energy, ρ
    else
        return total_energy
    end
end

function large_system_hartree_fock(mol::Molecule;
        max_iter = 100, δe = 1e-10,
        return_density = false
    )
    basis = Vector{STO3G}()
    Ne::Int = -mol.charge
    for atom::Atom in mol.atoms
        Ne += atom.Ne
        for orbital in atom.orbitals
            push!(basis, orbital)
        end
    end
    @assert iseven(Ne) "odd electron number is not suppport yet"
    Nb::Int = length(basis)
    Nn::Int = length(mol.atoms)
    
    # overlap matrix
    S = zeros(Float64, Nb, Nb) + I
    for i = 1:Nb
        for j = (i+1):Nb
            S[i, j] = S[j,i] = overlap_integral(basis[i], basis[j])
        end
    end
    
    # kinetic matrix
    T = zeros(Float64, Nb, Nb)
    for i = 1:Nb
        for j = i:Nb
            T[i,j] = T[j,i] = kinetic_integral(basis[i], basis[j])
        end
    end

    # nuclear elergy
    nuclear_energy = 0.0
    for I = 1:Nn
        for J = (I+1):Nn
            nuclear_energy += mol.atoms[I].Z * mol.atoms[J].Z / abs(mol.atoms[I].R - mol.atoms[J].R)
        end
    end

    # nuclear attraction matrix
    Vn = zeros(Float64, Nb, Nb)
    for atom in mol.atoms
        for i = 1:Nb
            for j = i:Nb
                _v = nuclear_acctaction_integral(atom.Z, atom.R, basis[i], basis[j])
                Vn[i,j] += _v
                if i != j
                    Vn[j, i] += _v
                end
            end
        end
    end

    # one body matrix
    H0 = T + Vn

    d, U = eigen(Symmetric(S))
    X = U * diagm(d.^(-1/2)) * U'

    ρ = zeros(Float64, Nb, Nb)  # density matrix
    G = zeros(Float64, Nb, Nb)  # effective interaction matrix    
    F = zeros(Float64, Nb, Nb)  # Fock matrix
    H = zeros(Float64, Nb, Nb)  # orthogonal Fock matrix
    C = zeros(Float64, Nb, Nb)  # orthogonal coefficients
    Cx = zeros(Float64, Nb, Nb) # real coefficients

    # prepare data for Veff! function, which links to libint2
    Nshell, l, alpha, coeff, R = make_shells(basis)

    electron_energy = 0.0
    total_energy = 0.0
    old_energy = Inf64

    # start hartree-fock iteration
    iter_times = 0
    for hf_iter = 1:max_iter
        Hx = H0 + G    
        H = X * Hx * X
        ε, C = eigen(H)
        Cx = real(X * C)

        # update density matrix
        ρ .= 0
        for i = 1:div(Ne, 2)
            for α = 1:Nb
                for a = 1:Nb
                    ρ[α, a] += 2*Cx[α, i] * Cx[a, i]
                end
            end
        end
        
        # update energy
        electron_energy = 0.0
        for α = 1:Nb
            for a = 1:Nb
                electron_energy += ρ[a, α] * (H0[α, a] + Hx[α, a])
            end
        end
        electron_energy *= 0.5

        if abs(old_energy - electron_energy) < δe
            break
        end

        old_energy = electron_energy

        # update effective interaction matrix
        Veff!(G, ρ, Nshell, l, alpha, coeff, R)
        iter_times = hf_iter
    end
    
    total_energy = electron_energy + nuclear_energy
    if iter_times == max_iter
        println("Warning: hartree_fock iterator time reaches `max_iter`, but may not reach the convergence condition")
    end
    if return_density
        return total_energy, ρ
    else
        return total_energy
    end
end