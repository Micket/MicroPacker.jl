using ArgParse

s = ArgParseSettings()
@add_arg_table! s begin
    "-f"
    help = "Output base filename"
    "-L"
    help = "Cell length (volume is L³)"
    arg_type = AbstractFloat
    default = 5.
    "-m"
    help = "Grid resolution. Total number of voxels is (m*L)³"
    arg_type = Int16
    default = Int16(10)
    "-s"
    help = "Seed for RNG for reproduceable runs"
    arg_type = Integer
    "--vol_frac_goal"
    help = "Goal for volume fraction inclusions (excluding overlap)"
    arg_type = Real
    default = 1.
    "--compress"
    help = "H5 compression (90 percent space savings, but not supported in DREAM.3D currently)"
    action = :store_true
end

add_arg_group!(s, "Packing")
@add_arg_table! s begin
    "--use_potential"
    help = "Use repulsive potential"
    action = :store_true
    "--nr_tries"
    help = "Number of random translations"
    arg_type = Integer
    default = 2500
    "--delta"
    help = "Maximum distance for randomized translations"
    default = 0.0
    arg_type = Real
    "--m_coarse"
    help = "Grid resolution during packing"
    arg_type = Int16
    default = Int16(10)
end

add_arg_group!(s, "WC grain shape")
@add_arg_table! s begin
    "-k"
    help = "k distribution"
    arg_type = String
    default = "U:0.4,1.4"
    "-r"
    help = "r distribution"
    arg_type = String
    default = "U:0.1,0.4"
    "-d"
    help = "d distribution"
    arg_type = String
    default = "U:0.5,1.5"
end

add_arg_group!(s, "Potts simulation")
@add_arg_table! s begin
    "--mc_steps"
    arg_type = AbstractFloat
    default = 0.05
    help = "Monte-Carlo steps (scales with (m*L)⁴. Set to zero to turn off."
    "--tau"
    help = "Ficticious temperature in Potts model"
    arg_type = AbstractFloat
    default = 0.5
end
