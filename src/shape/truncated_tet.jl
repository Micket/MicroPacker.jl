using StaticArrays

struct TruncatedTriangle{T}
    midpoint::SVector{3,T}
    rot::SMatrix{3,3,T,9}
    r::T
    k::T
    d_eq::T
    t::T
    vertices::SVector{6,SVector{2,T}}
end

function TruncatedTriangle(midpoint::SVector{3,T}, rot::SMatrix{3,3,T,9}, r::T, k::T, d_eq::T) where T
    L = cbrt(1/k * 8/3 * (2*r+1)^3 / ((r^2+4*r+1)*(r+1)) * 4*π/3) * (d_eq / 2)
    H = L * sqrt(T(3)) / 2
    h = (r+1)/(2*r+1) * H
    t = k * h

    a_long = 1/(2*r+1) * L
    a_short = r * a_long

    # vertices in coordinates of the triangle
    vertices = SA[
        SA[a_long/2, -H/3],
        SA[a_long/2 + a_short/2, sqrt(T(3))/2*a_short - H/3],
        SA[a_short/2, h - H/3],
        SA[-a_short/2, h - H/3],
        SA[-a_long/2 - a_short/2, sqrt(T(3))/2*a_short - H/3],
        SA[-a_long/2, -H/3]
    ]

    return TruncatedTriangle(midpoint, rot, r, k, d_eq, t, vertices)
end

struct TTGenerator{Fr,Fk,Fd}
    r::Fr
    k::Fk
    d_eq::Fd
end

function (g::TTGenerator)(midpoint::SVector{3,T}, rot::SMatrix{3,3,T,9}) where T
    TruncatedTriangle(midpoint, rot, g.r(), g.k(), g.d_eq())
end

function set_midpoint(t::TruncatedTriangle, midpoint::SVector{3,T}) where T
    return TruncatedTriangle(midpoint, t.rot, t.r, t.k, t.d_eq, t.t, t.vertices)
end

volume(shape::TruncatedTriangle) = π/6 * shape.d_eq^3

function bounding_box(shape::TruncatedTriangle)
    lower = shape.midpoint
    upper = shape.midpoint
    for vertex in shape.vertices
        vertex_a = shape.rot * SA[vertex[1], vertex[2], -shape.t/2]
        vertex_b = shape.rot * SA[vertex[1], vertex[2], shape.t/2]
        lower = min.(lower, vertex_a, vertex_b)
        upper = max.(upper, vertex_a, vertex_b)
    end
    return lower + shape.midpoint, upper + shape.midpoint
end

function inside(shape::TruncatedTriangle, p0::SVector{3,T}) where T
    testline(p, line) = (p[2] - line[1][2])*(line[2][1] - line[1][1]) - (line[2][2] - line[1][2])*(p[1] - line[1][1])

    # Rotate to coordinates of the triangle as used above
    r = shape.rot' * (p0 - shape.midpoint)

    # the triangle is within the x-y plane.
    return -shape.t/2 < r[3] < shape.t/2 &&
        testline(r, shape.vertices[SA[3,2]]) < 0 &&
        testline(r, shape.vertices[SA[4,3]]) < 0 &&
        testline(r, shape.vertices[SA[5,4]]) < 0 &&
        testline(r, shape.vertices[SA[5,6]]) > 0 &&
        testline(r, shape.vertices[SA[6,1]]) > 0 &&
        testline(r, shape.vertices[SA[1,2]]) > 0
#    return -shape.t/2 < r[3] < shape.t / 2 &&
#        testline(r, shape.vertices[SA[1,2]]) > 0 &&
#        testline(r, shape.vertices[SA[2,3]]) > 0 &&
#        testline(r, shape.vertices[SA[3,4]]) > 0 &&
#        testline(r, shape.vertices[SA[4,5]]) > 0 &&
#        testline(r, shape.vertices[SA[5,6]]) > 0 &&
#        testline(r, shape.vertices[SA[6,1]]) > 0 &&
end
