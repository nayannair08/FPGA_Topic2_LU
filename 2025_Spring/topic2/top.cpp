#include "dcl.h"

// Define a maximum number of neighbors per node.
#define MAX_NEIGHBORS 64

// ------------------------------
// Kernel Function: compute_tracks_HLS
// ------------------------------
// Reconstructs tracks from raw detector data using iterative DFS.
// "intt_required" indicates that only tracks with at least one hit with layer_id >= 3 are kept.
void compute_tracks_HLS(
    const int edge_index[2][MAX_EDGES],
    const data_t model_edge_probability[MAX_EDGES],
    const int num_edges,
    const int layer_id[MAX_HITS],
    const int n_pixels_arr[MAX_HITS],
    const data_t hit_cartesian[MAX_HITS][3],
    const int particle_id_arr[MAX_HITS],
    const data_t energy_arr[MAX_HITS],
    const data_t momentum_arr[MAX_HITS][3],
    const data_t track_origin_arr[MAX_HITS][3],
    const int trigger_node_arr[MAX_HITS],
    const int particle_type_arr[MAX_HITS],
    const int parent_particle_type_arr[MAX_HITS],
    const int num_hits,
    const data_t interaction_point_arr[3],
    bool trigger,
    bool has_trigger_pair,
    bool intt_required,
    EventInfo &event_info
) {
    // ----- Step 1: Filter Edges -----
    const data_t threshold = 0.5;
    int filtered_edges[2][MAX_EDGES];
    int filtered_count = 0;
    for (int i = 0; i < num_edges; i++) {
        if (model_edge_probability[i] >= threshold) {
            filtered_edges[0][filtered_count] = edge_index[0][i];
            filtered_edges[1][filtered_count] = edge_index[1][i];
            filtered_count++;
        }
    }
    
    // ----- Step 2: Build the Adjacency List -----
    int num_nodes = num_hits; // Each hit is a node.
    int adj_list[MAX_HITS][MAX_NEIGHBORS];
    int adj_count[MAX_HITS];
    // Initialize all neighbor counts to zero.
    for (int i = 0; i < num_nodes; i++) {
        adj_count[i] = 0;
    }
    // Fill the adjacency list using filtered edges.
    for (int i = 0; i < filtered_count; i++) {
        int u = filtered_edges[0][i];
        int v = filtered_edges[1][i];
        if (adj_count[u] < MAX_NEIGHBORS) {
            adj_list[u][adj_count[u]] = v;
            adj_count[u]++;
        }
        if (adj_count[v] < MAX_NEIGHBORS) {
            adj_list[v][adj_count[v]] = u;
            adj_count[v]++;
        }
    }
    
    // ----- Step 3: Connected Components via Iterative DFS -----
    bool visited[MAX_HITS];
    for (int i = 0; i < num_nodes; i++) {
        visited[i] = false;
    }

    int tracks[MAX_TRACKS][MAX_TRACK_SIZE];
    int track_sizes[MAX_TRACKS];
    int track_count = 0;

    // Iterate over every node. For each unvisited node, perform an iterative DFS.
    for (int i = 0; i < num_nodes; i++) {
        if (!visited[i]) {
            int current_track[MAX_TRACK_SIZE];
            int current_size = 0;
            int stack[MAX_HITS];
            int stack_ptr = 0;
            
            // Push the starting node.
            stack[stack_ptr++] = i;
            visited[i] = true;
            
            // DFS loop.
            while (stack_ptr > 0) {
                int node = stack[--stack_ptr];
                current_track[current_size++] = node;
                // Visit all neighbors.
                for (int j = 0; j < adj_count[node]; j++) {
                    int neigh = adj_list[node][j];
                    if (!visited[neigh]) {
                        visited[neigh] = true;
                        stack[stack_ptr++] = neigh;
                    }
                }
            }
            // Copy current_track into the global tracks array.
            for (int k = 0; k < current_size; k++) {
                tracks[track_count][k] = current_track[k];
            }
            track_sizes[track_count] = current_size;
            track_count++;
        }
    }
    
    // ----- Step 4: Process Each Track -----
    // Initialize EventInfo output buffers.
    for (int i = 0; i < MAX_TRACKS; i++) {
        for (int j = 0; j < NUM_LAYERS; j++) {
            event_info.n_pixels[i][j] = 0;
            event_info.track_n_hits[i][j] = 0;
            for (int k = 0; k < 3; k++) {
                event_info.track_hits[i][3*j + k] = 0;
            }
        }
        event_info.energy[i] = 0;
        for (int k = 0; k < 3; k++) {
            event_info.momentum[i][k] = 0;
            event_info.track_origin[i][k] = 0;
        }
        event_info.trigger_node[i] = 0;
        event_info.particle_id[i] = 0;
        event_info.particle_type[i] = 0;
        event_info.parent_particle_type[i] = 0;
    }
    
    int processed_tracks = 0;
    // For each track, optionally apply intt_required filtering.
    for (int t = 0; t < track_count && processed_tracks < MAX_TRACKS; t++) {
        if (intt_required) {
            bool valid_track = false;
            for (int k = 0; k < track_sizes[t]; k++) {
                int hit_idx = tracks[t][k];
                if (layer_id[hit_idx] >= 3) {
                    valid_track = true;
                    break;
                }
            }
            if (!valid_track) {
                continue; // Skip track if inner tracker requirement is not met.
            }
        }
        
        // Process per-layer group information.
        for (int j = 0; j < NUM_LAYERS; j++) {
            data_t weighted_sum[3] = {0, 0, 0};
            int sum_pixels = 0;
            int count_hits = 0;
            for (int k = 0; k < track_sizes[t]; k++) {
                int hit_idx = tracks[t][k];
                int hit_layer = layer_id[hit_idx];
                // Check if the hit falls within the current layer group.
                if (hit_layer >= layer_start[j] && hit_layer <= layer_end[j]) {
                    int pix = n_pixels_arr[hit_idx];
                    sum_pixels += pix;
                    for (int d = 0; d < 3; d++) {
                        weighted_sum[d] += pix * hit_cartesian[hit_idx][d];
                    }
                    count_hits++;
                }
            }
            for (int d = 0; d < 3; d++) {
                if (sum_pixels > 0)
                    event_info.track_hits[processed_tracks][3*j + d] = weighted_sum[d] / sum_pixels;
                else
                    event_info.track_hits[processed_tracks][3*j + d] = 0;
            }
            event_info.n_pixels[processed_tracks][j] = sum_pixels;
            event_info.track_n_hits[processed_tracks][j] = count_hits;
        }
        
        // Choose a representative hit (first hit) for overall track parameters.
        int rep_hit = tracks[t][0];
        event_info.energy[processed_tracks] = energy_arr[rep_hit];
        for (int d = 0; d < 3; d++) {
            event_info.momentum[processed_tracks][d] = momentum_arr[rep_hit][d];
            event_info.track_origin[processed_tracks][d] = track_origin_arr[rep_hit][d];
        }
        event_info.trigger_node[processed_tracks] = trigger_node_arr[rep_hit];
        event_info.particle_id[processed_tracks] = particle_id_arr[rep_hit];
        event_info.particle_type[processed_tracks] = particle_type_arr[rep_hit];
        event_info.parent_particle_type[processed_tracks] = parent_particle_type_arr[rep_hit];
        
        processed_tracks++;
    }
    event_info.num_tracks = processed_tracks;
    
    // ----- Step 5: Store Global Event Information -----
    for (int d = 0; d < 3; d++) {
        event_info.interaction_point[d] = interaction_point_arr[d];
    }
    event_info.trigger = trigger;
    event_info.has_trigger_pair = has_trigger_pair;
}
