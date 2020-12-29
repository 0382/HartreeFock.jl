# three dimension vector
import Base:+,-,*,/,^

struct Vec3{T<:Real}
    x::T
    y::T
    z::T
    Vec3{T}(x::T, y::T, z::T) where {T<:Real} = new(x, y, z)
end

Vec3(x::T, y::T, z::T) where {T<:Real} = Vec3{T}(x, y, z)
Vec3(x, y, z) = Vec3(promote(x, y, z)...)
vec3(x, y, z) = Vec3(x, y, z)
vec3_from_spherical(r::Float64, θ::Float64, φ::Float64) = Vec3(
    r*sin(θ)*cos(φ), r*sin(θ)*sin(φ), r*cos(θ)
)

Base.show(io::IO, v::Vec3) = begin
    print(io, "($(v.x), $(v.y), $(v.z))")
end

Base.eltype(::Vec3{T}) where {T<:Real} = T
Base.float(v::Vec3) = Vec3{Float64}(float(v.x), float(v.y), float(v.z))

+(v::Vec3) = v
-(v::Vec3) = Vec3(-v.x, -v.y, -v.z)
+(v1::Vec3, v2::Vec3) = Vec3(v1.x+v2.x, v1.y+v2.y, v1.z+v2.z)
-(v1::Vec3, v2::Vec3) = Vec3(v1.x-v2.x, v1.y-v2.y, v1.z-v2.z)
*(v1::Vec3, v2::Vec3) = v1.x*v2.x + v1.y*v2.y + v1.z*v2.z
*(k::Real, v::Vec3) = Vec3(k*v.x, k*v.y, k*v.z)
*(v::Vec3, k::Real) = k * v
/(v::Vec3, k::Real) = Vec3(v.x/k, v.y/k, v.z/k)
×(v1::Vec3, v2::Vec3) = Vec3(
    v1.y * v2.z - v1.z * v2.y,
    v1.z * v2.x - v1.x * v2.z,
    v1.x * v2.y - v1.y * v2.x
)

Base.abs2(v::Vec3) = v.x^2 + v.y^2 + v.z^2
Base.abs(v::Vec3) = hypot(v.x, v.y, v.z)
Base.length(v::Vec3) = abs(v)
distance(v1::Vec3, v2::Vec3) = abs(v1 - v2)

Base.getindex(v::Vec3, i::Integer) = begin
    i == 1 && return v.x
    i == 2 && return v.y
    i == 3 && return v.z
    error("out of range for Vec3 getindex i = $i")
end