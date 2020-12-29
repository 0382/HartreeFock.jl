function boys_function(n::Int, x::Float64)::Float64
    x < 1e-6 && return 1/(2n+1)
    n == 0 && return 0.5*√(π/x)*erf(√x)
    return ((2n-1)*boys_function(n-1, x) - exp(-x))/(2x)
end

"S(gs, gs)"
function overlap_integral(g1::Gaussian_s, g2::Gaussian_s)::Float64
    p = (g1.α + g2.α)
    μ = (g1.α * g2.α) / p
    γ = 4μ / p
    return γ^(3/4) * exp(-μ * abs2(g1.R - g2.R))
end

"S(gs, gp)"
function overlap_integral(g1::Gaussian_s, g2::Gaussian_p)::Float64
    p = g1.α + g2.α
    μ = (g1.α * g2.α) / p
    γ = 4μ / p
    Rab = g1.R - g2.R
    return √(g1.α) * γ^(5/4) * exp(-μ * abs2(Rab)) * (Rab * g2.n)
end

"S(gp, gs)"
function overlap_integral(g1::Gaussian_p, g2::Gaussian_s)::Float64
    overlap_integral(g2, g1)
end

"S(gp, gp)"
function overlap_integral(g1::Gaussian_p, g2::Gaussian_p)::Float64
    p = g1.α + g2.α
    μ = (g1.α * g2.α) / p
    γ = 4μ / p
    Rab = g1.R - g2.R
    return γ^(5/4) * ((g1.n * g2.n) -√(g1.α*g2.α*γ) * (Rab * g1.n) * (Rab * g2.n)) * exp(-μ * abs2(Rab))
end

"T(gs, gs)"
function kinetic_integral(g1::Gaussian_s, g2::Gaussian_s)::Float64
    p = g1.α + g2.α
    μ = (g1.α * g2.α) / p
    γ = 4μ / p
    Rab2 = abs2(g1.R - g2.R)
    return γ^(3/4) * μ * (3-2μ*Rab2) * exp(-μ * Rab2)
end

"T(gp, gs)"
function kinetic_integral(g1::Gaussian_p, g2::Gaussian_s)::Float64
    p = g1.α + g2.α
    μ = g1.α * g2.α / p
    γ = 4μ / p
    Rab = g1.R - g2.R
    Rab2 = abs2(Rab)
    return γ^(5/4) * exp(-μ * Rab2) * μ * (2μ*Rab2-5) * √g2.α * (Rab * g1.n)
end

"T(gs, gp)"
function kinetic_integral(g1::Gaussian_s, g2::Gaussian_p)::Float64
    kinetic_integral(g2, g1)
end

"T(gp, gp)"
function kinetic_integral(g1::Gaussian_p, g2::Gaussian_p)::Float64
    p = g1.α + g2.α
    μ = g1.α * g2.α / p
    γ = 4μ / p
    Rab = g1.R - g2.R
    dotnab = g1.n * g2.n
    μRab2 = μ * abs2(Rab)
    M = (g1.n * Rab) * (g2.n * Rab)
    return γ^(5/4) * exp(-μRab2) * μ * (5dotnab - 14μ * M - 2μRab2 * dotnab + 4μ * μRab2 * M)
end

"Vn(gs, gs)"
function nuclear_acctaction_integral(Zc::Int, Rc::Vec3, g1::Gaussian_s, g2::Gaussian_s)::Float64
    p = g1.α + g2.α
    μ = (g1.α * g2.α) / p
    γ = 4μ / p
    Rμ = (g1.α / p) * g1.R + (g2.α / p) * g2.R
    return -2Zc * γ^(3/4) * √(p/π) * exp(-μ * abs2(g1.R - g2.R)) * boys_function(0, p*abs2(Rμ-Rc))
end

"Vn(gp,gs)"
function nuclear_acctaction_integral(Zc::Int, Rc::Vec3, g1::Gaussian_p, g2::Gaussian_s)::Float64
    p = g1.α + g2.α
    μ = (g1.α * g2.α) / p
    γ = 4μ / p
    Rμ = (g1.α / p) * g1.R + (g2.α / p) * g2.R
    pR2 = p*abs2(Rμ - Rc)
    return -4Zc * γ^(3/4) * √(g1.α * p / π) * exp(-μ * abs2(g1.R-g2.R)) *
        (-(g1.n * (Rμ-Rc)) * boys_function(1, pR2) + (g1.n * (Rμ-g1.R)) * boys_function(0, pR2))
end

"Vn(gs,gp)"
function nuclear_acctaction_integral(Zc::Int, Rc::Vec3, g1::Gaussian_s, g2::Gaussian_p)::Float64
    nuclear_acctaction_integral(Zc, Rc, g2, g1)
end

"Vn(gp, gp)"
function nuclear_acctaction_integral(Zc::Int, Rc::Vec3, g1::Gaussian_p, g2::Gaussian_p)::Float64
    p = g1.α + g2.α
    μ = g1.α * g2.α / p
    γ = 4μ / p
    Rab = g1.R - g2.R
    Rμ = (g1.α / p) * g1.R + (g2.α / p) * g2.R
    Rμc = Rμ - Rc
    Dac = g1.n * Rμc
    Dbc = g2.n * Rμc
    Daa = g1.n * (Rμ - g1.R)
    Dbb = g2.n * (Rμ - g2.R)
    x = p * abs2(Rμc)
    F0x = boys_function(0, x)
    F1x = boys_function(1, x)
    F2x = boys_function(2, x)
    return -2Zc * γ^(5/4) * √(p/π) * exp(-μ*abs2(g1.R - g2.R)) * (
        (g1.n * g2.n) * (F0x - F1x) + 2p * Dac * Dbc * F2x -
        2p * (Dac * Dbb + Dbc * Daa) * F1x + 2p * Daa * Dbb * F0x
    )
end

function electron_repulsion_integral(g1::Gaussian_s, g2::Gaussian_s, g3::Gaussian_s, g4::Gaussian_s)
    p1 = g1.α + g4.α
    p2 = g2.α + g3.α
    μ1 = (g1.α * g4.α) / p1
    μ2 = (g2.α * g3.α) / p2
    γ1 = 4μ1 / p1
    γ2 = 4μ2 / p2
    Rμ1 = (g1.α / p1) * g1.R + (g4.α / p1) * g4.R
    Rμ2 = (g2.α / p2) * g2.R + (g3.α / p2) * g3.R
    λ = p1 * p2 / (p1 + p2)
    R = abs(Rμ1 - Rμ2)
    if R < 1e-6
        return (γ1*γ2)^(3/4) * exp(-μ1 * abs2(g1.R-g4.R)-μ2 * abs2(g2.R-g3.R)) * 2 * √(λ / π)
    else
        return (γ1*γ2)^(3/4) * exp(-μ1 * abs2(g1.R-g4.R)-μ2 * abs2(g2.R-g3.R)) * erf(√λ * R) / R
    end
end

function two_center_integral(s1::STONG{N}, s2::STONG{N}, integral::Function) where N
    total = 0.0
    for i = 1:N
        for j = 1:N
            total += s1.d[i] * s2.d[j] * integral(s1.g[i], s2.g[j])
        end
    end
    return total
end

function four_center_integral(s1::STONG{N}, s2::STONG{N}, s3::STONG{N}, s4::STONG{N}, integral::Function) where N
    total = 0.0
    for i = 1:N
        for j = 1:N
            for k = 1:N
                for l = 1:N
                    total += s1.d[i]*s2.d[j]*s3.d[k]*s4.d[l] * integral(s1.g[i], s2.g[j], s3.g[k], s4.g[l])
                end
            end
        end
    end
    return total
end

overlap_integral(s1::STONG, s2::STONG) = two_center_integral(s1, s2, overlap_integral)
kinetic_integral(s1::STONG, s2::STONG) = two_center_integral(s1, s2, kinetic_integral)
nuclear_acctaction_integral(Zc::Int, Rc::Vec3, s1::STONG, s2::STONG) = two_center_integral(
    s1, s2,
    (g1,g2)->nuclear_acctaction_integral(Zc, Rc, g1, g2)
)
electron_repulsion_integral(s1::STONG, s2::STONG, s3::STONG, s4::STONG) = four_center_integral(
    s1, s2, s3, s4,
    electron_repulsion_integral
)