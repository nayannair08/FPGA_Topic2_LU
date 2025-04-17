open_project project_1

# set top function of the HLS design
set_top compute_tracks_HLS

# add source file
add_files top.cpp

# add testbench
add_files -tb host.cpp

# add data file
add_files -tb event100607116/edge_index.bin
add_files -tb event100607116/edge_phi_slope.bin
add_files -tb event100607116/edge_z0.bin
add_files -tb event100607116/energy.bin
add_files -tb event100607116/has_trigger_pair.bin
add_files -tb event100607116/hit_cartesian.bin
add_files -tb event100607116/hit_cylindrical.bin
add_files -tb event100607116/interaction_point.bin
add_files -tb event100607116/layer_id.bin
add_files -tb event100607116/model_edge_probability.bin
add_files -tb event100607116/momentum.bin
add_files -tb event100607116/n_pixels.bin
add_files -tb event100607116/parent_particle_type.bin
add_files -tb event100607116/particle_id.bin
add_files -tb event100607116/particle_type.bin
add_files -tb event100607116/phi_slope_max.bin
add_files -tb event100607116/track_origin.bin
add_files -tb event100607116/trigger.bin
add_files -tb event100607116/trigger_node.bin
add_files -tb event100607116/z0_max.bin

open_solution "solution8"

# FPGA part and clock configuration
set_part {xczu3eg-sbva484-1-e}

# default frequency is 100 MHz
#create_clock -period 4 -name default

# C synthesis for HLS design, generating RTL
#csynth_design

# C/RTL co-simulation; can be commented if not needed
#cosim_design

# export generated RTL as an IP; can be commented if not needed
export_design -format ip_catalog -flow impl



exit