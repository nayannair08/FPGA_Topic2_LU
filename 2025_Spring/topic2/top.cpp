#include "dcl.h"

// Define a maximum number of neighbors per node.
#define MAX_NEIGHBORS 64

// -----------------------
// Recursive DFS Function
// -----------------------
// This function recursively visits nodes in the graph, building a connected component.
// Parameters:
//   node         : current node to process.
//   current_track: array to store nodes of the current connected component (track).
//   track_size   : reference to the current track size (number of nodes recorded so far).
//   visited      : boolean array marking which nodes have been visited.
//   adj_list     : adjacency list for each node.
//   adj_count    : number of neighbors for each node.
void dfs_recursive(
    int node,
    int current_track[MAX_TRACK_SIZE],
    int &track_size,
    bool visited[MAX_HITS],
    int adj_list[MAX_HITS][MAX_NEIGHBORS],
    int adj_count[MAX_HITS]
) {
    visited[node] = true;
    current_track[track_size++] = node;
    for (int j = 0; j < adj_count[node]; j++) {
        int neigh = adj_list[node][j];
        if (!visited[neigh]) {
            dfs_recursive(neigh, current_track, track_size, visited, adj_list, adj_count);
        }
    }
}

// ------------------------------
// Kernel Function: compute_tracks_HLS
// ------------------------------
// Reconstructs tracks from raw detector data using recursive DFS.
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
    int num_nodes = num_hits; // Each hit corresponds to a node.
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
    
    // ----- Step 3: Connected Components via Recursive DFS -----
    bool visited[MAX_HITS];
    for (int i = 0; i < num_nodes; i++) {
        visited[i] = false;
    }
    
    int tracks[MAX_TRACKS][MAX_TRACK_SIZE];
    int track_sizes[MAX_TRACKS];
    int track_count = 0;
    
    // For each node, if it hasn't been visited, start a recursive DFS.
    for (int i = 0; i < num_nodes; i++) {
        if (!visited[i]) {
            int current_track[MAX_TRACK_SIZE];
            int current_size = 0;
            dfs_recursive(i, current_track, current_size, visited, adj_list, adj_count);
            // Copy current_track to the global tracks array.
            for (int k = 0; k < current_size; k++) {
                tracks[track_count][k] = current_track[k];
            }
            track_sizes[track_count] = current_size;
            track_count++;
        }
    }
    
    // ----- Step 4: Process Each Track -----
    // Zero out the output EventInfo buffers.
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
    // For each reconstructed track, compute weighted averages and assign track parameters.
    for (int t = 0; t < track_count && t < MAX_TRACKS; t++) {
        // For each detector layer group.
        for (int j = 0; j < NUM_LAYERS; j++) {
            data_t weighted_sum[3] = {0, 0, 0};
            int sum_pixels = 0;
            int count_hits = 0;
            // Loop over the hits in the current track.
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
            // Compute weighted average hit position for this layer group.
            for (int d = 0; d < 3; d++) {
                if (sum_pixels > 0)
                    event_info.track_hits[processed_tracks][3*j + d] = weighted_sum[d] / sum_pixels;
                else
                    event_info.track_hits[processed_tracks][3*j + d] = 0;
            }
            event_info.n_pixels[processed_tracks][j] = sum_pixels;
            event_info.track_n_hits[processed_tracks][j] = count_hits;
        }
        
        // Choose a representative hit for overall track parameters.
        // Here, we simply choose the first hit in the track.
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
