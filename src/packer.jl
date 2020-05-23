using StaticArrays
using LinearAlgebra
using Optim

function random_axis(t::Type{T}) where T <: AbstractFloat
    θ = rand(t)*2*π
    ϕ = acos(2*rand(t)-1)
    return SA[sin(ϕ)*cos(θ), sin(ϕ)*sin(θ), cos(ϕ)]
end

function random_rotation(t::Type{T}) where T <: AbstractFloat
    small = T(0.1)
    r0 = random_axis(t)
    dum = zero(T)
    while dum < small
        r1 = random_axis(t)
        r1 -= (r1 ⋅ r0)*r0
        dum = norm(r1)
    end
    r1 /= dum
    r2 = r0 × r1
    return hcat(r0, r1, r2)
end

function euler_angles(R)
    cos_θ = R[3, 3]
    sin_θ_sq = 1 - cos_θ^2

    if sin_θ_sq < 1e-11
        θ = acos(cos_θ)
        ψ = zero(eltype(R))
        ϕ = mod(atan(R[2, 1], R[1, 1]), 2π)
        return SA[ϕ, θ, ψ]
    end

    θ = acos(cos_θ)
    ϕ = mod(atan(R[1, 3], -R[2, 3]), 2π)
    ψ = mod(atan(R[3, 1], R[3, 2]), 2π)

    return SA[ϕ, θ, ψ]
end

function prepare_inclusions(vol_frac_goal::AbstractFloat, L::AbstractFloat, inclusiontype::Type, generator)
    ftype = typeof(L)
    total_volume = L^3

    # Make random grains
    total_inclusion_volume = zero(ftype)
    inclusions = inclusiontype[]
    while total_inclusion_volume < vol_frac_goal * total_volume
        midpoint = L * @SVector rand(ftype, 3)
        rot = random_rotation(ftype)
        inclusion = generator(midpoint, rot)
        push!(inclusions, inclusion)
        total_inclusion_volume += volume(inclusion)
    end

    return inclusions
end

function voxelize(shape, L::AbstractFloat, M::Integer)
    delta_x = L/M

    # Inexpensive bounding box of prism to limit the voxels to check:
    lower, upper = bounding_box(shape)
    loweri = floor.(typeof(M), lower*M/L)
    upperi = ceil.(typeof(M), upper*M/L)

    voxel_indices = Array{SVector{3, typeof(M)},1}(undef, prod(Int32.(upperi - loweri)))
    total_voxels = 0
    for iz in loweri[3]:upperi[3], iy in loweri[2]:upperi[2], ix in loweri[1]:upperi[1]
        i = SA[ix, iy, iz]
        r = delta_x .* (0.5 .+ i)
        if inside(shape, r)
            total_voxels += 1
            voxel_indices[total_voxels] = i
        end
    end
    return voxel_indices[1:total_voxels]
end

function pack_inclusions!(inclusions::AbstractArray, M::T, L::AbstractFloat, nr_tries::Integer, Δ::T) where T<:Integer
    grain_count_type = Int16
    grain_ids = zeros(grain_count_type, M, M, M)
    overlaps = zeros(grain_count_type, M, M, M)
    inclusions_voxels = Array{SVector{3, T},1}[]
    inc_voxels = 0
    M3 = Int64(M)^3
    for (i, inclusion) in enumerate(inclusions)
        inclusion_voxels = voxelize(inclusion, L, M)

        overlap_min = M3
        Δ_min = zeros(SVector{3, typeof(M)})
        Δ_try = Δ_min
        for n_try in 1:nr_tries
            overlap_try = 0
            for inclusion_voxel in inclusion_voxels
                wrapped_voxel = one(typeof(M)) .+ mod.(inclusion_voxel + Δ_try, M)
                if grain_ids[wrapped_voxel...] > 1 # claimed, so add overlap
                    overlap_try += 1
                end
            end

            if overlap_try < overlap_min
                overlap_min = overlap_try
                Δ_min = Δ_try
                if overlap_min == 0
                    break
                end
            end
            Δ_try = rand(-Δ:Δ, SVector{3})
        end

        # Move the truncated triangle to the right position as well:
        inclusions[i] = set_midpoint(inclusion, inclusion.midpoint + mod.(Δ_min * L/M, L))

        # Rerun with optimal delta_x, delta_y, delta_z
        for (j, inclusion_voxel) in enumerate(inclusion_voxels)
            wrapped_voxel = one(typeof(M)) .+ mod.(inclusion_voxel + Δ_min, M)
            inclusion_voxels[j] = wrapped_voxel
            if grain_ids[wrapped_voxel...] == 1 # still unclaimed binder
                grain_ids[wrapped_voxel...] = i+2
                inc_voxels += 1
            elseif grain_ids[wrapped_voxel...] > 1 # claimed, so add overlap
                overlaps[wrapped_voxel...] += 1
            end
        end
        push!(inclusions_voxels, sort(inclusion_voxels))
        println("grain $i: WC fraction: $(inc_voxels), tries: $nr_tries delta: $Δ_min")
    end

    return grain_ids, overlaps, inclusions_voxels
end

phases(grain_ids) = grain_ids .== 0

# Computes euler angles as used by Dream3D
function euler_angles(grain_ids, inclusions)
    inclusions_angles = euler_angles.(trunc_triangles)

    euler_angles = zeros(size(grain_ids)..., 3)
    for i in CartesianIndices(grain_ids)
        if grain_ids[i] > 1
            grain_id = grain_ids[i] - 1
            euler_angles[i, :] = inclusions_angles[i]
        end
    end

    return euler_angles
end

# Helper
function check_nb(id, old_id, new_id)
    if id == old_id
        return 1 # the new id is different, add 1 area
    elseif id == new_id
        return -1 # the new id is the same, subtract 1 area
    else
        return 0
    end
end

function neighbours(i::CartesianIndex{3}, M::Integer)
    # right, left, forward, backward, up, down
    sub(x) = ifelse(x > 1, x-1, M)
    add(x) = ifelse(x < M, x+1, 1)
    return SA[
        CartesianIndex(sub(i.I[1]), i.I[2], i.I[3]),
        CartesianIndex(add(i.I[1]), i.I[2], i.I[3]),
        CartesianIndex(i.I[1], sub(i.I[2]), i.I[3]),
        CartesianIndex(i.I[1], add(i.I[2]), i.I[3]),
        CartesianIndex(i.I[1], i.I[2], sub(i.I[3])),
        CartesianIndex(i.I[1], i.I[2], add(i.I[3])),
    ]
end

# Monte-Carlo Potts model suitable for minimizing ground boundaries by simulation grain boundary migrations.
# Effectively minimizes contiguity
function make_mcp_bound!(grain_ids, gb_voxels, overlaps, voxel_indices, steps::Integer, kBT::Real)
    M = size(grain_ids)[1] # TODO allow non-square domains

    # Compute the exponentials for different number of neighbouring voxels
    exp_ΔA_kBT = exp.(-SA[1,2,3,4]/kBT)

    overlap_index = findall((x)-> x > 0, overlaps)
    if len(overlap_index) == 0
        return
    end

    for step in 1:steps
        # Choose a random voxel in the overlapping regions
        gb_voxel_index = rand(overlap_index)

        # Check if it is a gb voxel
        if gb_voxels[gb_voxel_index]
            gb_voxel_id = grain_ids[gb_voxel_index]
            nb_indices = neighbours(gb_voxel_index, M)

            # Clear and populate the *set* of different neighbor indices that are allowed; id > 1, i.e. not binder
            nb_ids = grain_ids[nb_indices]

            # Find all possible values the voxel could switch to.
            nb_set = @SVector zeros(6)
            nr_diff_ids = 0
            for i in 1:6
                nb_id = nb_ids[i]
                if nb_id > 1 && nb_id != gb_voxel_id && nb_id ∈ nb_set
                    # Not in set and gb_voxel_id belongs to nb_id and can thus be changed to nb_id
                    if insorted(voxel_indices[nb_id-1], gb_voxel_index)
                        nr_diff_ids += 1
                        nb_set[nr_diff_ids] = nb_id
                    end
                end
            end

            # the set of allowed changes can be zero if we are at the edges of a truncated triangle
            if nr_diff_ids > 0
                new_id = nb_set[rand(1:nr_diff_ids)]

                ΔA = check_nb.(nb_ids, gb_voxel_id, new_id)
                sum_ΔA = sum(ΔA)

                # Metropolis algorithm
                if sum_ΔA <= 0 || (sum_ΔA > 0 && rand() < exp_ΔA_kBT[sum_ΔA-1])
                    grain_ids[gb_voxel_index] = new_id
                    gb_voxels[gb_voxel_index] += sum_ΔA
                    for (δ, nb_index) in zip(ΔA, nb_indices)
                        gb_voxels[nb_index] += δ
                    end
                end
            end
        end
    end
end
