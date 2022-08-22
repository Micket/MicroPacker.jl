using StaticArrays

struct Cylinder{T}
    midpoint::SVector{3,T}
    rot::SMatrix{3,3,T,9}
    r::T
    h::T
end

struct CylinderGenerator{Fr,Fh}
    r::Fr
    h::Fh
end

function (g::CylinderGenerator)(midpoint::SVector{3,T}, rot::SMatrix{3,3,T,9}) where T
    Cylinder(midpoint, rot, g.r(), g.h())
end

function set_midpoint(shape::Cylinder, midpoint::SVector{3,T}) where T
    return Cylinder(midpoint, shape.rot, shape.r, shape.h)
end

volume(shape::Cylinder) = π * shape.r^2 * h

function bounding_box(shape::Cylinder)
    lower = (shape.midpoint + shape.rot * SA[0, 0, -shape.h/2]) .- shape.r
    upper = (shape.midpoint + shape.rot * SA[0, 0,  shape.h/2]) .+ shape.r
    return lower, upper
end

function inside(shape::Cylinder, p0::SVector{3,T}) where T
    p = shape.rot' * (p0 - shape.midpoint)
    return -shape.h/2 < p[3] < shape.h/2 && shape.r^2 < (p[1]^2 + p[2]^2)
end
