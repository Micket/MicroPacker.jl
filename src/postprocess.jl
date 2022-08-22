using StaticArrays
using LinearAlgebra
using Optim

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

mass_fraction(vol_frac, ρ_i, ρ_m) = vol_frac*ρ_i / (vol_frac*ρ_i + (1-vol_frac)*ρ_m)
