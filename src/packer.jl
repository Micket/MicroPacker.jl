using StaticArrays
using LinearAlgebra

function random_axis(t::Type{T}) where T
    θ = rand(t)*2*π
    ϕ = acos(2*rand(t)-1)
    return SA[sin(ϕ)*cos(θ), sin(ϕ)*sin(θ), cos(ϕ)]
end

function random_rotation(t::Type{T}) where T
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

function prepare_inclusions(vol_frac_goal::Real, L::Real, inclusiontype::Type, generator)
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

function voxelize(shape, L::Real, M::Integer)
    delta_x = L/M

    # Inexpensive bounding box to limit the voxels to check:
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

function pack_inclusions!(inclusions::AbstractArray, M::T, L::Real, nr_tries::Integer, Δ::T) where T<:Integer
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
