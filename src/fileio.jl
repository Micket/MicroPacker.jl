using HDF5

function write_dream3d(filename, L, grain_ids, phases, good_voxels, euler_angles;
               surface_voxels=nothing, gb_voxels=nothing, interface_voxels=nothing, overlaps=nothing, compress=false)
    M = size(grain_ids)
    spacing = L ./ size(grain_ids)
    h5open(filename + ".dream3d", 'w') do f
        cmpr = compress ? "gzip" : nothing

        attrs(f)["FileVersion"] = "7.0"
        attrs(f)["DREAM3D Version"] = "6.2"

        grp_containers = g_create(f, "DataContainers")
        grp_voxel_data = g_create(grp_containers, "VoxelDataContainer")

        attrs(grp_voxel_data)["NUM_POINTS"] = prod(M)
        attrs(grp_voxel_data)["VTK_DATA_OBJECT"] = "VTK_STRUCTURED_POINTS"

        # Geometry group
        grp_simpl_geometry = g_create(grp_voxel_data, "_SIMPL_GEOMETRY")
        attrs(grp_simpl_geometry)["GeometryName"] = "ImageGeometry"
        attrs(grp_simpl_geometry)["GeometryTypeName"] = "ImageGeometry"
        attrs(grp_simpl_geometry)["GeometryType"] = 11
        attrs(grp_simpl_geometry)["SpatialDimensionality"] = 3
        attrs(grp_simpl_geometry)["UnitDimensionality"] = 3
        dset_dimensions = d_create(grp_simpl_geometry, "DIMENSIONS", datatype(Int64), (3,))
        dset_dimensions[:] = M

        dset_origin = d_create(grp_simpl_geometry, "ORIGIN", datatype(Float32), (3,))
        dset_origin[:] = [0, 0, 0]

        dset_origin = d_create(grp_simpl_geometry, "SPACING", datatype(Float32), (3,))
        dset_origin[:] = spacing

        # Cell data group
        grp_cell_data = g_create(grp_voxel_data, "CellData")
        attrs(grp_cell_data)["AttributeMatrixType"] = [3]
        attrs(grp_cell_data)["TupleDimensions"] = M

        dset_grain_ids = d_create(grp_cell_data, "GrainIds", datatype(Int32), (M[0]*M[1]*M[2],), compression=cmpr)
        attrs(dset_grain_ids)["ComponentDimensions"] = [1]
        attrs(dset_grain_ids)["DataArrayVersion"] = [2]
        attrs(dset_grain_ids)["ObjectType"] = "DataArray<int32_t>"
        attrs(dset_grain_ids)["TupleDimensions"] = M
        dset_grain_ids[:] = grain_ids

        dset_phases = d_create(grp_cell_data, "Phases", datatype(Int32), (M[0]*M[1]*M[2],), compression=cmpr)
        attrs(dset_phases)["ComponentDimensions"] = [1]
        attrs(dset_phases)["DataArrayVersion"] = [2]
        attrs(dset_phases)["ObjectType"] = "DataArray<int32_t>"
        attrs(dset_phases)["TupleDimensions"] = M
        dset_phases[:] = phases

        dset_good_voxels = d_create(grp_cell_data, "GoodVoxels", datatype(UInt8), (M[0]*M[1]*M[2],), compression=cmpr)
        attrs(dset_good_voxels)["ComponentDimensions"] = [1]
        attrs(dset_good_voxels)["DataArrayVersion"] = [2]
        attrs(dset_good_voxels)["ObjectType"] = "DataArray<bool>"
        attrs(dset_good_voxels)["TupleDimensions"] = M
        dset_good_voxels[:] = good_voxels

        dset_euler_angles = d_create(grp_cell_data, "EulerAngles", datatype(Float32), (M[0]*M[1]*M[2], 3), compression=cmpr)
        attrs(dset_euler_angles)["ComponentDimensions"] = [3]
        attrs(dset_euler_angles)["DataArrayVersion"] = [2]
        attrs(dset_euler_angles)["ObjectType"] = "DataArray<float>"
        attrs(dset_euler_angles)["TupleDimensions"] = M
        dset_euler_angles[:] = euler_angles

        # Cell ensemble data
        grp_ensemble_data = g_create(grp_voxel_data, "CellEnsembleData")
        attrs(grp_ensemble_data)["AttributeMatrixType"] = [11]
        attrs(grp_ensemble_data)["Name"] = "EnsembleData"
        attrs(grp_ensemble_data)["TupleDimensions"] = [2]

        dset_crystal_structs = d_create(grp_ensemble_data, "CrystalStructures", datatype(UInt32), (3,))
        attrs(dset_crystal_structs)["ComponentDimensions"] = [1]
        attrs(dset_crystal_structs)["DataArrayVersion"] = [2]
        attrs(dset_crystal_structs)["ObjectType"] = "DataArray<uint32_t>"
        attrs(dset_crystal_structs)["TupleDimensions"] = [2]
        dset_crystal_structs[:] = [999, 1, 0]

        dset_phase_types = d_create(grp_ensemble_data, "PhaseTypes", datatype(UInt32), (3,))
        attrs(dset_phase_types)["ComponentDimensions"] = [1]
        attrs(dset_phase_types)["DataArrayVersion"] = [2]
        attrs(dset_phase_types)["ObjectType"] = "DataArray<uint32_t>"
        attrs(dset_phase_types)["TupleDimensions"] = [2]
        dset_phase_types[:] = [999, 3, 1]
        if surface_voxels != nothing
            dset_surface_voxels = d_create(grp_cell_data, "SurfaceVoxels", datatype(Int8), (M[0]*M[1]*M[2],), compression=cmpr)
            attrs(dset_surface_voxels)["DataArrayVersion"] = [2]
            attrs(dset_surface_voxels)["TupleDimensions"] = M
            attrs(dset_surface_voxels)["ComponentDimensions"] = [1]
            attrs(dset_surface_voxels)["ObjectType"] = "DataArray<int8_t>"
            dset_surface_voxels[:] = surface_voxels
        end
        if gb_voxels != nothing
            dset_gb_voxels = d_create(grp_cell_data, "GrainBoundaryVoxels", datatype(Int8), (M[0]*M[1]*M[2],), compression=cmpr)
            attrs(dset_gb_voxels)["DataArrayVersion"] = [2]
            attrs(dset_gb_voxels)["TupleDimensions"] = M
            attrs(dset_gb_voxels)["ComponentDimensions"] = [1]
            attrs(dset_gb_voxels)["ObjectType"] = "DataArray<int8_t>"
            dset_gb_voxels[:] = gb_voxels
        end
        if interface_voxels != nothing
            dset_interface_voxels = d_create(grp_cell_data, "InterfaceVoxels", datatype(Int8), (M[0]*M[1]*M[2],), compression=cmpr)
            attrs(dset_interface_voxels)["DataArrayVersion"] = [2]
            attrs(dset_interface_voxels)["TupleDimensions"] = M
            attrs(dset_interface_voxels)["ComponentDimensions"] = [1]
            attrs(dset_interface_voxels)["ObjectType"] = "DataArray<int8_t>"
            dset_interface_voxels[:] = interface_voxels
        end
        if overlaps != nothing
            dset_overlaps = d_create(grp_cell_data, "Overlaps", datatype(Int8), (M[0]*M[1]*M[2],), compression=cmpr)
            attrs(dset_overlaps)["DataArrayVersion"] = [2]
            attrs(dset_overlaps)["TupleDimensions"] = M
            attrs(dset_overlaps)["ComponentDimensions"] = [1]
            attrs(dset_overlaps)["ObjectType"] = "DataArray<int8_t>"
            dset_overlaps[:] = overlaps
        end

        pipeline = g_create(f, "Pipeline")
        attrs(pipeline)["Number_Filters"] = [0]
    end

    write_xdmf(filename, M, spacing, surface_voxels != Nothing, gb_voxels != Nothing,
               interface_voxels != Nothing, overlaps != Nothing)
end

function write_xdmf(filename, M, spacing; surface_voxels=false, gb_voxels=false, interface_voxels=false, overlaps=false)
    open(filename + ".xdmf", 'w') do f
        filename_hdf5 = rsplit(filename, '/')[end] + ".dream3d"
        dimensions = "$(M[1]) $(M[2]) $(M[3])"

        write(f, "<?xml version=\"1.0\"?>\n")
        write(f, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\"[]>\n")
        write(f, "<Xdmf xmlns:xi=\"http://www.w3.org/2003/XInclude\" Version=\"2.2\">\n")
        write(f, " <Domain>\n\n")
        write(f, "  <Grid Name=\"Cell Data\" GridType=\"Uniform\">\n")
        write(f, "    <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"$(M[1]+1) $(M[2]+1) $(M[3]+1) \"></Topology>\n")
        write(f, "    <Geometry Type=\"ORIGIN_DXDYDZ\">\n")
        write(f, "      <!-- Origin -->\n")
        write(f, "      <DataItem Format=\"XML\" Dimensions=\"3\">0 0 0</DataItem>\n")
        write(f, "      <!-- DxDyDz (Spacing/Resolution)-->\n")
        write(f, "      <DataItem Format=\"XML\" Dimensions=\"3\">$(spacing[1]) $(spacing[2]) $(spacing[3])</DataItem>\n")
        write(f, "    </Geometry>\n")
        write(f, "    <Attribute Name=\"EulerAngles (Cell)\" AttributeType=\"Vector\" Center=\"Cell\">\n")
        write(f, "      <DataItem Format=\"HDF\" Dimensions=\"$(dimensions) 3\" NumberType=\"Float\" Precision=\"4\" >\n")
        write(f, "        $(filename_hdf5):/DataContainers/VoxelDataContainer/CellData/EulerAngles\n")
        write(f, "      </DataItem>\n")
        write(f, "    </Attribute>\n")
        write(f, "\n")
        write(f, "    <Attribute Name=\"GoodVoxels (Cell)\" AttributeType=\"Scalar\" Center=\"Cell\">\n")
        write(f, "      <DataItem Format=\"HDF\" Dimensions=\"$(dimensions)\" NumberType=\"uchar\" Precision=\"1\" >\n")
        write(f, "        $(filename_hdf5):/DataContainers/VoxelDataContainer/CellData/GoodVoxels\n")
        write(f, "      </DataItem>\n")
        write(f, "    </Attribute>\n")
        write(f, "\n")
        write(f, "    <Attribute Name=\"GrainIds (Cell)\" AttributeType=\"Scalar\" Center=\"Cell\">\n")
        write(f, "      <DataItem Format=\"HDF\" Dimensions=\"$(dimensions)\" NumberType=\"Int\" Precision=\"4\" >\n")
        write(f, "        $(filename_hdf5):/DataContainers/VoxelDataContainer/CellData/GrainIds\n")
        write(f, "      </DataItem>\n")
        write(f, "    </Attribute>\n")
        write(f, "\n")
        write(f, "    <Attribute Name=\"Phases (Cell)\" AttributeType=\"Scalar\" Center=\"Cell\">\n")
        write(f, "      <DataItem Format=\"HDF\" Dimensions=\"$(dimensions)\" NumberType=\"Int\" Precision=\"4\" >\n")
        write(f, "        $(filename_hdf5):/DataContainers/VoxelDataContainer/CellData/Phases\n")
        write(f, "      </DataItem>\n")
        write(f, "    </Attribute>\n")
        write(f, "\n")
        if surface_voxels
            write(f, "    <Attribute Name=\"SurfaceVoxels (Cell)\" AttributeType=\"Scalar\" Center=\"Cell\">\n")
            write(f, "      <DataItem Format=\"HDF\" Dimensions=\"$(dimensions)\" NumberType=\"Char\" Precision=\"1\" >\n")
            write(f, "        $(filename_hdf5):/DataContainers/VoxelDataContainer/CellData/SurfaceVoxels\n")
            write(f, "      </DataItem>\n")
            write(f, "    </Attribute>\n")
            write(f, "\n")
        end
        if gb_voxels
            write(f, "    <Attribute Name=\"GrainBoundaryVoxels (Cell)\" AttributeType=\"Scalar\" Center=\"Cell\">\n")
            write(f, "      <DataItem Format=\"HDF\" Dimensions=\"$(dimensions)\" NumberType=\"Char\" Precision=\"1\" >\n")
            write(f, "        $(filename_hdf5):/DataContainers/VoxelDataContainer/CellData/GrainBoundaryVoxels\n")
            write(f, "      </DataItem>\n")
            write(f, "    </Attribute>\n")
            write(f, "\n")
        end
        if interface_voxels
            write(f, "    <Attribute Name=\"InterfaceVoxels (Cell)\" AttributeType=\"Scalar\" Center=\"Cell\">\n")
            write(f, "      <DataItem Format=\"HDF\" Dimensions=\"$(dimensions)\" NumberType=\"Char\" Precision=\"1\" >\n")
            write(f, "        $(filename_hdf5):/DataContainers/VoxelDataContainer/CellData/InterfaceVoxels\n")
            write(f, "      </DataItem>\n")
            write(f, "    </Attribute>\n")
            write(f, "\n")
        end
        if overlaps
            write(f, "    <Attribute Name=\"Overlaps (Cell)\" AttributeType=\"Scalar\" Center=\"Cell\">\n")
            write(f, "      <DataItem Format=\"HDF\" Dimensions=\"$(dimensions)\" NumberType=\"Char\" Precision=\"1\" >\n")
            write(f, "        $(filename_hdf5):/DataContainers/VoxelDataContainer/CellData/Overlaps\n")
            write(f, "      </DataItem>\n")
            write(f, "    </Attribute>\n")
            write(f, "\n")
        end

        write(f, "  </Grid>\n")
        write(f, "    <!-- *************** END OF Cell Data *************** -->\n")
        write(f, "\n")
        write(f, " </Domain>\n")
        write(f, "</Xdmf>\n")
    end
end

function read_dream3d(filename)
    h5open(filename + ".dream3d", 'r') do f
        # To return:
        # M, spacing, grain_ids, phases, good_voxels, euler_angles, overlaps=None
        grp_voxel_data = f["DataContainers"]["VoxelDataContainer"]
        grp_simpl_geometry = grp_voxel_data["_SIMPL_GEOMETRY"]
        M = read(grp_simpl_geometry, "DIMENSIONS")
        spacing = read(grp_simpl_geometry, "SPACING")

        grp_cell_data = grp_voxel_data["CellData"]
        grain_ids = read(grp_cell_data, "GrainIds")
        phases = read(grp_cell_data, "Phases")
        good_voxels = read(grp_cell_data, "GoodVoxels")
        euler_angles = read(grp_cell_data, "EulerAngles")
        if exists(grp_cell_data, "Overlaps")
            overlaps = read(grp_cell_data, "Overlaps")
        end
    end

    return M, spacing, grain_ids, phases, good_voxels, euler_angles, overlaps
end
