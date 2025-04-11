#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cassert>

// Define data structures
std::vector<std::vector<int> > edge_index;
std::vector<float> model_edge_probability;
std::vector<float> energy;
std::vector<float> hit_cylindrical;
std::vector<int> parent_particle_type;
std::vector<float> phi_slope_max;
std::vector<int> trigger_node;
std::vector<float> edge_phi_slope;
std::vector<bool> has_trigger_pair;
std::vector<float> interaction_point;
std::vector<float> momentum;
std::vector<int> particle_id;
std::vector<float> track_origin;
std::vector<float> z0_max;
std::vector<float> edge_z0;
std::vector<float> hit_cartesian;
std::vector<int> layer_id;
std::vector<int> n_pixels;
std::vector<int> particle_type;
std::vector<bool> trigger;

// General dimensions (to be derived from file size)
int num_hits, num_edges;

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

// ---------- Logging Functions ----------
template<typename T>
void log_vector(const std::string& name, const std::vector<T>& vec, int max_items = 10) {
    std::cout << name << " (size=" << vec.size() << "): ";
    for (size_t i = 0; i < vec.size() && i < max_items; ++i) {
        std::cout << vec[i] << " ";
    }
    if (vec.size() > max_items) std::cout << "...";
    std::cout << std::endl;
}

template<typename T>
void log_vector_2d(const std::string& name, const std::vector<std::vector<T> >& vec2d, int max_rows = 2, int max_cols = 10) {
    std::cout << name << " (shape=" << vec2d.size() << "x" << (vec2d.empty() ? 0 : vec2d[0].size()) << "):" << std::endl;
    for (size_t i = 0; i < vec2d.size() && i < max_rows; ++i) {
        std::cout << "  [";
        for (size_t j = 0; j < vec2d[i].size() && j < max_cols; ++j) {
            std::cout << vec2d[i][j] << (j + 1 < vec2d[i].size() ? ", " : "");
        }
        if (vec2d[i].size() > max_cols) std::cout << "...";
        std::cout << "]" << std::endl;
    }
    if (vec2d.size() > max_rows) std::cout << "  ..." << std::endl;
}

void log_dataset() {
    std::cout << "\n--- Logging Dataset ---\n";
    log_vector_2d("edge_index", edge_index);
    log_vector("model_edge_probability", model_edge_probability);
    log_vector("energy", energy);
    log_vector("hit_cylindrical", hit_cylindrical);
    log_vector("parent_particle_type", parent_particle_type);
    log_vector("phi_slope_max", phi_slope_max);
    log_vector("trigger_node", trigger_node);
    log_vector("edge_phi_slope", edge_phi_slope);
    log_vector("has_trigger_pair", has_trigger_pair);
    log_vector("interaction_point", interaction_point);
    log_vector("momentum", momentum);
    log_vector("particle_id", particle_id);
    log_vector("track_origin", track_origin);
    log_vector("z0_max", z0_max);
    log_vector("edge_z0", edge_z0);
    log_vector("hit_cartesian", hit_cartesian);
    log_vector("layer_id", layer_id);
    log_vector("n_pixels", n_pixels);
    log_vector("particle_type", particle_type);
    log_vector("trigger", trigger);
    std::cout << "--- End Logging ---\n";
}

// ---------- Main ----------
int main() {
    std::string directory_path = "/Users/sahithsr/Desktop/Gatech/Spring_2025/FPG/Project/FPGA_Topic2_LU/2025_Spring/topic2/event000041188/"; // Set your path here

    // Infer number of hits and edges
    num_hits = infer_size<float>(directory_path + "hit_cartesian.bin", 3);
    num_edges = infer_size<int>(directory_path + "edge_index.bin", 2);

    // Read edge_index as 2D
    edge_index.resize(2, std::vector<int>(num_edges));
    {
        std::ifstream file(directory_path + "edge_index.bin", std::ios::binary);
        assert(file.is_open());
        for (int i = 0; i < 2; ++i) {
            file.read(reinterpret_cast<char*>(edge_index[i].data()), num_edges * sizeof(int));
        }
        file.close();
    }

    read_binary_file(directory_path + "model_edge_probability.bin", model_edge_probability);
    read_binary_file(directory_path + "energy.bin", energy);
    read_binary_file(directory_path + "hit_cylindrical.bin", hit_cylindrical);
    read_binary_file(directory_path + "parent_particle_type.bin", parent_particle_type);
    read_binary_file(directory_path + "phi_slope_max.bin", phi_slope_max);
    read_binary_file(directory_path + "trigger_node.bin", trigger_node);
    read_binary_file(directory_path + "edge_phi_slope.bin", edge_phi_slope);

    std::vector<uint8_t> has_trigger_pair_raw;
    read_binary_file(directory_path + "has_trigger_pair.bin", has_trigger_pair_raw);
    has_trigger_pair.resize(has_trigger_pair_raw.size());
    for (size_t i = 0; i < has_trigger_pair_raw.size(); ++i)
        has_trigger_pair[i] = static_cast<bool>(has_trigger_pair_raw[i]);

    read_binary_file(directory_path + "interaction_point.bin", interaction_point);
    read_binary_file(directory_path + "momentum.bin", momentum);
    read_binary_file(directory_path + "particle_id.bin", particle_id);
    read_binary_file(directory_path + "track_origin.bin", track_origin);
    read_binary_file(directory_path + "z0_max.bin", z0_max);
    read_binary_file(directory_path + "edge_z0.bin", edge_z0);
    read_binary_file(directory_path + "hit_cartesian.bin", hit_cartesian);
    read_binary_file(directory_path + "layer_id.bin", layer_id);
    read_binary_file(directory_path + "n_pixels.bin", n_pixels);
    read_binary_file(directory_path + "particle_type.bin", particle_type);

    std::vector<uint8_t> trigger_raw;
    read_binary_file(directory_path + "trigger.bin", trigger_raw);
    trigger.resize(trigger_raw.size());
    for (size_t i = 0; i < trigger_raw.size(); ++i)
        trigger[i] = static_cast<bool>(trigger_raw[i]);

    std::cout << "Binary files read successfully." << std::endl;
    std::cout << "Number of hits: " << num_hits << std::endl;
    std::cout << "Number of edges: " << num_edges << std::endl;

    log_dataset(); // 👈 Log everything

    return 0;
}

