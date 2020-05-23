include("utils.jl")
#include("fileio.jl")
include("packer.jl")
include("shape/truncated_tet.jl")
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
sort!(inclusions, by=volume, rev=true)

println("Prepared $(length(inclusions)) inclusions")

#tt = TruncatedTriangle(SA[0., 0., 0.], SA[1. 0. 0.; 0. 1. 0.; 0. 0. 1.], r₀, k₀, d₀)
#lower, upper = bounding_box(tt)
#test = voxelize(inclusions[1], L, M)

grain_ids, overlaps, inclusions_voxels = pack_inclusions!(inclusions, M, L, nr_tries, Δ)

m_phase = phases(grain_ids)
m_euler_angles = euler_angles(grain_ids, inclusions)

make_mcp_bound!(grain_ids, inclusions_voxels, overlaps, voxel_indices, steps::Integer, kBT::Real)

write_dream3d(options["f"], SA[L, L, L], grain_ids, phases, euler_angles)
