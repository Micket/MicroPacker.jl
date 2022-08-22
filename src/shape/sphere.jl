using StaticArrays

struct Sphere{T}
    midpoint::SVector{3,T}
    r::T
end

struct SphereGenerator{Fr}
    r::Fr
end

function (g::SphereGenerator)(midpoint::SVector{3,T}, rot) where T
    Sphere(midpoint, g.r())
end

function set_midpoint(shape::Sphere, midpoint::SVector{3,T}) where T
    return Sphere(midpoint, shape.r)
end

volume(shape::Sphere) = 4/3 * π * shape.r^3

function bounding_box(shape::Sphere)
    lower = shape.midpoint .- shape.r
    upper = shape.midpoint .+ shape.r
    return lower, upper
end

function inside(shape::Sphere, p0::SVector{3,T}) where T
    p = p0 - shape.midpoint
    return shape.r^2 < sum(p.^2)
end
