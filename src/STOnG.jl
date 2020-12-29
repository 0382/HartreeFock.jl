# STO-nG basis definition

abstract type BasisFunction end
abstract type GTO <: BasisFunction end
abstract type GaussianFunction <: GTO end

"Gaussian s wave function"
struct Gaussian_s <: GaussianFunction
    α::Float64       # gussian function exponent
    R::Vec3{Float64} # center (or actually nuclues' position)
end

"Gaussian p wave function"
struct Gaussian_p <:GaussianFunction
    α::Float64       # gussian function exponent
    R::Vec3{Float64} # center
    n::Vec3{Float64} # orbital direction
end

"STO-nG function"
struct STONG{N} <: GTO
    d::NTuple{N, Float64}           # coefficients
    g::NTuple{N, GaussianFunction}  # gaussian functions
    STONG{N}(d::NTuple{N, Float64}, g::NTuple{N, GaussianFunction}) where N = new(d, g)
end

STONG(d::NTuple{N, Float64}, g::NTuple{N, GaussianFunction}) where N = STONG{N}(d, g)

struct STONG_coeff{N}
    coefficients::NTuple{N, Float64}
    exponents::NTuple{N,Float64}
    STONG_coeff{N}(coeff::NTuple{N, Float64}, _exp::NTuple{N, Float64}) where N = new(coeff, _exp)
end

"construct STONG from STONG_coeff"
STONG(coeff::STONG_coeff{N}, center::Vec3{Float64}; type::Symbol=:s) where N = begin
    if type === :s
        return STONG{N}(
            coeff.coefficients,
            ntuple(i->Gaussian_s(coeff.exponents[i], center), N)
        )
    elseif type === :px
        return STONG{N}(
            coeff.coefficients,
            ntuple(i->Gaussian_p(coeff.exponents[i], center, vec3(1., 0., 0.)), N)
        )
    elseif type === :py
        return STONG{N}(
            coeff.coefficients,
            ntuple(i->Gaussian_p(coeff.exponents[i], center, vec3(0., 1., 0.)), N)
        )
    elseif type === :pz
        return STONG{N}(
            coeff.coefficients,
            ntuple(i->Gaussian_p(coeff.exponents[i], center, vec3(0., 0., 1.)), N)
        )
    end
    error("orbital type $type is not support yet!")
end

const STO3G = STONG{3}
const STO3G_coeff = STONG_coeff{3}

# the data comes from basis_set_exchange
# https://github.com/MolSSI-BSE/basis_set_exchange
const coefficient_data = Dict(
    "H" => Dict(
        "1s" => STO3G_coeff(
            (0.1543289673E+00, 0.5353281423E+00, 0.4446345422E+00), # coefficients
            (0.3425250914E+01, 0.6239137298E+00, 0.1688554040E+00)  # exponents
        )
    ),
    "He" => Dict(
        "1s" => STO3G_coeff(
            (0.1543289673E+00, 0.5353281423E+00, 0.4446345422E+00), # coefficients
            (0.6362421394E+01, 0.1158922999E+01, 0.3136497915E+00)  # exponents
        )
    ),
    "C" => Dict(
        "1s" => STO3G_coeff(
            (0.1543289673E+00, 0.5353281423E+00, 0.4446345422E+00), # coefficients
            (0.7161683735E+02, 0.1304509632E+02, 0.3530512160E+01)  # exponents
        ),
        "2s" => STO3G_coeff(
            (-0.9996722919E-01, 0.3995128261E+00, 0.7001154689E+00),# coefficients
            (0.2941249355E+01, 0.6834830964E+00, 0.2222899159E+00)  # exponents
        ),
        "2p" => STO3G_coeff(
            (0.1559162750E+00, 0.6076837186E+00, 0.3919573931E+00), # coefficients
            (0.2941249355E+01, 0.6834830964E+00, 0.2222899159E+00)  # exponents
        )
    ),
    "N" => Dict(
        "1s" => STO3G_coeff(
            (0.1543289673E+00, 0.5353281423E+00, 0.4446345422E+00), # coefficients
            (0.9910616896E+02, 0.1805231239E+02, 0.4885660238E+01)  # exponents
        ),
        "2s" => STO3G_coeff(
            (-0.9996722919E-01, 0.3995128261E+00, 0.7001154689E+00),# coefficients
            (0.3780455879E+01, 0.8784966449E+00, 0.2857143744E+00)  # exponents
        ),
        "2p" => STO3G_coeff(
            (0.1559162750E+00, 0.6076837186E+00, 0.3919573931E+00), # coefficients
            (0.3780455879E+01, 0.8784966449E+00, 0.2857143744E+00)  # exponents
        )
    ),
    "O" => Dict(
        "1s" => STO3G_coeff(
            (0.1543289673E+00, 0.5353281423E+00, 0.4446345422E+00), # coefficients
            (0.1307093214E+03, 0.2380886605E+02, 0.6443608313E+01)  # exponents
        ),
        "2s" => STO3G_coeff(
            (-0.9996722919E-01, 0.3995128261E+00, 0.7001154689E+00),# coefficients
            (0.5033151319E+01, 0.1169596125E+01, 0.3803889600E+00)  # exponents
        ),
        "2p" => STO3G_coeff(
            (0.1559162750E+00, 0.6076837186E+00, 0.3919573931E+00), # coefficients
            (0.5033151319E+01, 0.1169596125E+01, 0.3803889600E+00)  # exponents
        )
    ),
    "F" => Dict(
        "1s" => STO3G_coeff(
            (0.1543289673E+00, 0.5353281423E+00, 0.4446345422E+00), # coefficients
            (0.1666791340E+03, 0.3036081233E+02, 0.8216820672E+01)  # exponents
        ),
        "2s" => STO3G_coeff(
            (-0.9996722919E-01, 0.3995128261E+00, 0.7001154689E+00),# coefficients
            (0.6464803249E+01, 0.1502281245E+01, 0.4885884864E+00)  # exponents
        ),
        "2p" => STO3G_coeff(
            (0.1559162750E+00, 0.6076837186E+00, 0.3919573931E+00), # coefficients
            (0.6464803249E+01, 0.1502281245E+01, 0.4885884864E+00)  # exponents
        )
    )
)

Base.length(::STONG{N}) where N = N

Base.show(io::IO, gs::Gaussian_s) = begin
    print(io, "gs(α=$(gs.α), R=$(gs.R))")
end

Base.show(io::IO, gp::Gaussian_p) = begin
    print(io, "gp(α=$(gp.α), R=$(gp.R), n=$(gp.n))")
end

Base.show(io::IO, stong::STONG{N}) where N = begin
    println(io, "STO-$(N)G")
    for i = 1:N
        println(io, "  $(stong.d[i]) × $(stong.g[i])")
    end
end