include("utils.jl")
#include("fileio.jl")
include("packer.jl")
include("truncated_tet.jl")
include("args.jl")

options = parse_args(ARGS, s)

if options["s"] != nothing
    Random.seed!(options['s'])
end

d_eq, d₀ = parse_dist(options["d"])
r, r₀ = parse_dist(options["r"])
k, k₀ = parse_dist(options["k"])
generator = TTGenerator(r, k, d_eq)

vol_frac_goal = options["vol_frac_goal"]
L = 2. #options["L"]
m = 5 #options["m"]
M = round(Int16, m * L / d₀)
nr_tries = 5 #options["nr_tries"]
Δ = M÷Int8(2) #options["delta"]
spacing = L / M

inclusions = prepare_inclusions(vol_frac_goal, L, TruncatedTriangle, generator)

# Sort triangles w.r.t. volume typically leads to better packing:
sort!(inclusions, by=x->x.d_eq, rev=true)

println("Prepared $(length(inclusions)) inclusions")

if options["use_potential"]
    println("use_potential: TODO")
    #ccb.optimize_midpoints(L, trunc_triangles)
end

#tt = TruncatedTriangle(SA[0., 0., 0.], SA[1. 0. 0.; 0. 1. 0.; 0. 0. 1.], r₀, k₀, d₀)
#lower, upper = bounding_box(tt)
#test = voxelize(inclusions[1], L, M)

grain_ids, overlaps, inclusions_voxels = pack_inclusions!(inclusions, M, L, nr_tries, Δ)

phases = grain_ids .== 0
good_voxels = zeros(Bool, size(grain_ids))
euler_angles = zeros((size(grain_ids)..., 3))

make_mcp_bound!(grain_ids, inclusions_voxels, overlaps, voxel_indices, steps::Integer, kBT::Real)

write_dream3d(options["f"], SA[L, L, L], grain_ids, phases, good_voxels, euler_angles)
