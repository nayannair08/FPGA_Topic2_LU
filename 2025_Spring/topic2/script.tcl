open_project project_3

# set top function of the HLS design
set_top compute_tracks_HLS

# add source file
add_files top.cpp 

# add testbench
add_files -tb host.cpp 

# add data file
# grab every .bin in any event directory two levels under ./converted
set all_event_bins [glob -nocomplain ./event_dir/*/event*/*.bin]

# add them as test-bench data files
add_files -tb $all_event_bins

open_solution "solution2"

#config_compile -cflags "-std=c++17" -lflags "-lstdc++fs"

# FPGA part and clock configuration
set_part {xczu3eg-sbva484-1-e}

# default frequency is 100 MHz
#create_clock -period 4 -name default

# C synthesis for HLS design, generating RTL
csynth_design

# C/RTL co-simulation; can be commented if not needed
cosim_design

# export generated RTL as an IP; can be commented if not needed
export_design -format ip_catalog -flow impl



exit