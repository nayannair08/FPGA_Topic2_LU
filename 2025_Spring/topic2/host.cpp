#include "dcl.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cassert>

// ---------- Binary File Reader ----------
template<typename T>
void read_binary_file(const std::string& filename, std::vector<T>& container) {
    std::ifstream file(filename, std::ios::binary);
    assert(file.is_open());
    file.seekg(0, std::ios::end);
    size_t size = file.tellg();
    file.seekg(0);
    container.resize(size / sizeof(T));
    file.read(reinterpret_cast<char*>(container.data()), size);
    file.close();
}

template<typename T>
int infer_size(const std::string& filename, int columns = 1) {
    std::ifstream file(filename, std::ios::binary);
    assert(file.is_open());
    file.seekg(0, std::ios::end);
    size_t size = file.tellg();
    file.close();
    return size / (sizeof(T) * columns);
}

int main() {
    // Set the directory path where the binary files reside.
    //std::string directory_path = "/nethome/sbommu3/FPGA/project/FPGA_Topic2_LU/2025_Spring/topic2/event000041188/"; // Adjust to your directory
    std::string directory_path = "/nethome/sbommu3/FPGA/project/FPGA_Topic2_LU/2025_Spring/topic2/event100607116/"; // Adjust to your directory
    // Raw data containers (vectors to temporarily hold file data)
    std::vector<int> edge_index_raw;
    std::vector<data_t> model_edge_prob;
    std::vector<int> layer_id_vec;
    std::vector<int> n_pixels_vec;
    std::vector<data_t> hit_cartesian_vec;
    std::vector<int> particle_id_vec;
    std::vector<data_t> energy_vec;
    std::vector<data_t> momentum_vec;
    std::vector<data_t> track_origin_vec;
    std::vector<int> trigger_node_vec;
    std::vector<int> particle_type_vec;
    std::vector<int> parent_particle_type_vec;
    std::vector<data_t> interaction_point_vec;
    std::vector<int> trigger_vec;  // Assuming stored as int: 0 or 1
    std::vector<uint8_t> has_trigger_pair_raw;

    // Determine number of edges and hits using file sizes
    int num_edges = infer_size<int>(directory_path + "edge_index.bin", 2);
    int num_hits  = infer_size<data_t>(directory_path + "hit_cartesian.bin", 3);

    std::cout << "Number of edges: " << num_edges << std::endl;
    std::cout << "Number of hits: " << num_hits << std::endl;

    // Read binary files
    read_binary_file(directory_path + "edge_index.bin", edge_index_raw);
    read_binary_file(directory_path + "model_edge_probability.bin", model_edge_prob);
    read_binary_file(directory_path + "layer_id.bin", layer_id_vec);
    read_binary_file(directory_path + "n_pixels.bin", n_pixels_vec);
    read_binary_file(directory_path + "hit_cartesian.bin", hit_cartesian_vec);
    read_binary_file(directory_path + "particle_id.bin", particle_id_vec);
    read_binary_file(directory_path + "energy.bin", energy_vec);
    read_binary_file(directory_path + "momentum.bin", momentum_vec);
    read_binary_file(directory_path + "track_origin.bin", track_origin_vec);
    read_binary_file(directory_path + "trigger_node.bin", trigger_node_vec);
    read_binary_file(directory_path + "particle_type.bin", particle_type_vec);
    read_binary_file(directory_path + "parent_particle_type.bin", parent_particle_type_vec);
    read_binary_file(directory_path + "interaction_point.bin", interaction_point_vec);
    read_binary_file(directory_path + "trigger.bin", trigger_vec);
    read_binary_file(directory_path + "has_trigger_pair.bin", has_trigger_pair_raw);

    // Convert trigger flag and has_trigger_pair to bool (nonzero is true)
    bool trigger_flag = (trigger_vec.size() > 0 && trigger_vec[0] != 0);
    bool has_trigger_pair = (has_trigger_pair_raw.size() > 0 && has_trigger_pair_raw[0] != 0);
    
    // Define intt_required flag (set to true if an inner tracker hit is required).
    bool intt_required = true;

    // ---------- Pack Data into Fixed-Size Arrays ----------
    int edge_index_arr[2][MAX_EDGES];
    for (int i = 0; i < num_edges; i++) {
        edge_index_arr[0][i] = edge_index_raw[2 * i];
        edge_index_arr[1][i] = edge_index_raw[2 * i + 1];
    }

    data_t model_edge_probability_arr[MAX_EDGES];
    for (int i = 0; i < num_edges; i++) {
        model_edge_probability_arr[i] = model_edge_prob[i];
    }

    int layer_id_arr[MAX_HITS];
    int n_pixels_arr[MAX_HITS];
    data_t hit_cartesian_arr[MAX_HITS][3];
    int particle_id_arr[MAX_HITS];
    data_t energy_arr[MAX_HITS];
    data_t momentum_arr[MAX_HITS][3];
    data_t track_origin_arr[MAX_HITS][3];
    int trigger_node_arr[MAX_HITS];
    int particle_type_arr[MAX_HITS];
    int parent_particle_type_arr[MAX_HITS];

    for (int i = 0; i < num_hits; i++) {
        layer_id_arr[i] = layer_id_vec[i];
        n_pixels_arr[i] = n_pixels_vec[i];
    }
    for (int i = 0; i < num_hits; i++) {
        for (int j = 0; j < 3; j++) {
            hit_cartesian_arr[i][j] = hit_cartesian_vec[3 * i + j];
        }
    }
    for (int i = 0; i < num_hits; i++) {
        particle_id_arr[i] = particle_id_vec[i];
    }
    for (int i = 0; i < num_hits; i++) {
        energy_arr[i] = energy_vec[i];
    }
    for (int i = 0; i < num_hits; i++) {
        for (int j = 0; j < 3; j++) {
            momentum_arr[i][j] = momentum_vec[3 * i + j];
        }
    }
    for (int i = 0; i < num_hits; i++) {
        for (int j = 0; j < 3; j++) {
            track_origin_arr[i][j] = track_origin_vec[3 * i + j];
        }
    }
    for (int i = 0; i < num_hits; i++) {
        trigger_node_arr[i] = trigger_node_vec[i];
    }
    for (int i = 0; i < num_hits; i++) {
        particle_type_arr[i] = particle_type_vec[i];
        parent_particle_type_arr[i] = parent_particle_type_vec[i];
    }
    data_t interaction_point_arr[3];
    for (int i = 0; i < 3; i++) {
        interaction_point_arr[i] = interaction_point_vec[i];
    }

    // ---------- Call the Track Reconstruction Kernel ----------
    EventInfo event_info; // This structure will be filled by the kernel.
    compute_tracks_HLS(
        edge_index_arr,
        model_edge_probability_arr,
        num_edges,
        layer_id_arr,
        n_pixels_arr,
        hit_cartesian_arr,
        particle_id_arr,
        energy_arr,
        momentum_arr,
        track_origin_arr,
        trigger_node_arr,
        particle_type_arr,
        parent_particle_type_arr,
        num_hits,
        interaction_point_arr,
        trigger_flag,
        has_trigger_pair,
        intt_required,  // Pass the new inner tracker requirement flag.
        event_info
    );

    // ---------- Output Results ----------
    std::cout << "\n--- Reconstructed Event Information ---" << std::endl;
    // Global event-level information.
    std::cout << "Global Interaction Point: [" << event_info.interaction_point[0] << ", "
              << event_info.interaction_point[1] << ", " 
              << event_info.interaction_point[2] << "]" << std::endl;
    std::cout << "Trigger: " << event_info.trigger << std::endl;
    std::cout << "Has Trigger Pair: " << event_info.has_trigger_pair << std::endl;
    std::cout << "Number of Tracks: " << event_info.num_tracks << std::endl;

    // Print details for each track.
    for (int i = 0; i < event_info.num_tracks; i++) {
        std::cout << "\nTrack " << i << ":" << std::endl;
        std::cout << "  Energy: " << event_info.energy[i] << std::endl;
        std::cout << "  Momentum: [" 
                  << event_info.momentum[i][0] << ", " 
                  << event_info.momentum[i][1] << ", " 
                  << event_info.momentum[i][2] << "]" << std::endl;
        std::cout << "  Track Origin: [" 
                  << event_info.track_origin[i][0] << ", " 
                  << event_info.track_origin[i][1] << ", " 
                  << event_info.track_origin[i][2] << "]" << std::endl;
        std::cout << "  Trigger Node: " << event_info.trigger_node[i] << std::endl;
        std::cout << "  Particle ID: " << event_info.particle_id[i] << std::endl;
        std::cout << "  Particle Type: " << event_info.particle_type[i] << std::endl;
        std::cout << "  Parent Particle Type: " << event_info.parent_particle_type[i] << std::endl;
        // Print per-layer-group information.
        for (int j = 0; j < NUM_LAYERS; j++) {
            std::cout << "  Layer Group " << j << ":" << std::endl;
            std::cout << "    n_pixels: " << event_info.n_pixels[i][j] << std::endl;
            std::cout << "    track_n_hits: " << event_info.track_n_hits[i][j] << std::endl;
            std::cout << "    track_hits: ["
                      << event_info.track_hits[i][3*j] << ", "
                      << event_info.track_hits[i][3*j+1] << ", "
                      << event_info.track_hits[i][3*j+2] << "]" << std::endl;
        }
    }
    return 0;
}
