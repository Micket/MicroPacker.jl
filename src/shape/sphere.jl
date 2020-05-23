using StaticArrays

struct Sphere{T}
    midpoint::SVector{3,T}
    r::T
end

struct SphereGenerator{Fr}
    r::Fr
end

function (g::SphereGenerator)(midpoint, rot)
    Sphere(midpoint, g.r())
end

volume(shape::Sphere) = 4/3 * π * shape.r^3

function set_midpoint(shape::Sphere, midpoint)
    return Sphere(midpoint, t.r)
end

function bounding_box(shape::Sphere)
    lower = shape.midpoint .- shape.r
    upper = shape.midpoint .+ shape.r
    return lower, upper
end

function inside(shape::Sphere, p0)
    p = p0 - shape.midpoint
    return shape.r^2 < sum(p.^2)
end
