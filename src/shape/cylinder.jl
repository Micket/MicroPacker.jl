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

function (g::CylinderGenerator)(midpoint, rot)
    Cylinder(midpoint, rot, g.r(), g.h())
end

volume(shape::Cylinder) = π * shape.r^2 * h

function set_midpoint(shape::Cylinder, midpoint)
    return Cylinder(midpoint, t.rot, t.r, t.h)
end

function bounding_box(shape::Cylinder)
    lower = (shape.midpoint + shape.rot * SA[0, 0, -shape.h/2]) .- shape.r
    upper = (shape.midpoint + shape.rot * SA[0, 0,  shape.h/2]) .+ shape.r
    return lower, upper
end

function inside(shape::Cylinder, p0)
    p = shape.rot' * (p0 - shape.midpoint)
    return -shape.h/2 < p[3] < shape.h/2 && shape.r^2 < (p[1]^2 + p[2]^2)
end
