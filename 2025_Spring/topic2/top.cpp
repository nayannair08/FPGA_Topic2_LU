#include "dcl.h"
#include <map>
#include <cmath>
#define DEBUG

#define MAX_NEIGHBORS 64

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

    int num_nodes = num_hits;
    int adj_list[MAX_HITS][MAX_NEIGHBORS];
    int adj_count[MAX_HITS];
    for (int i = 0; i < num_nodes; i++) adj_count[i] = 0;
    for (int i = 0; i < filtered_count; i++) {
        int u = filtered_edges[0][i];
        int v = filtered_edges[1][i];
        if (adj_count[u] < MAX_NEIGHBORS) adj_list[u][adj_count[u]++] = v;
        if (adj_count[v] < MAX_NEIGHBORS) adj_list[v][adj_count[v]++] = u;
    }

    bool visited[MAX_HITS] = {false};
    int tracks[MAX_TRACKS][MAX_TRACK_SIZE];
    int track_sizes[MAX_TRACKS];
    int track_count = 0;

    for (int i = 0; i < num_nodes; i++) {
        if (!visited[i]) {
            int current_track[MAX_TRACK_SIZE], current_size = 0;
            int stack[MAX_HITS], stack_ptr = 0;
            stack[stack_ptr++] = i;
            visited[i] = true;

            while (stack_ptr > 0) {
                int node = stack[--stack_ptr];
                current_track[current_size++] = node;
                for (int j = 0; j < adj_count[node]; j++) {
                    int neigh = adj_list[node][j];
                    if (!visited[neigh]) {
                        visited[neigh] = true;
                        stack[stack_ptr++] = neigh;
                    }
                }
            }
            for (int k = 0; k < current_size; k++) {
                tracks[track_count][k] = current_track[k];
            }
            track_sizes[track_count++] = current_size;
        }
    }

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
            if (!valid_track) continue;
        }

        for (int j = 0; j < NUM_LAYERS; j++) {
            data_t weighted_sum[3] = {0};
            int sum_pixels = 0, count_hits = 0;
            for (int k = 0; k < track_sizes[t]; k++) {
                int hit_idx = tracks[t][k];
                int hit_layer = layer_id[hit_idx];
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
                event_info.track_hits[processed_tracks][3*j + d] = sum_pixels > 0 ? weighted_sum[d] / sum_pixels : 0;
            }
            event_info.n_pixels[processed_tracks][j] = sum_pixels;
            event_info.track_n_hits[processed_tracks][j] = count_hits;
        }

        std::map<int, int> pid_counts;
        int max_count = 0;
        int selected_hit_idx = tracks[t][0];

        for (int k = 0; k < track_sizes[t]; k++) {
            int idx = tracks[t][k];
            int pid = particle_id_arr[idx];
            if (std::isnan(static_cast<double>(pid))) continue;
            pid_counts[pid]++;
            if (pid_counts[pid] > max_count) {
                max_count = pid_counts[pid];
                selected_hit_idx = idx;
            }
        }

        if (max_count == 0) {
            for (int k = 0; k < track_sizes[t]; k++) {
                int idx = tracks[t][k];
                if (std::isnan(static_cast<double>(particle_id_arr[idx]))) {
                    selected_hit_idx = idx;
                    break;
                }
            }
        }

        event_info.energy[processed_tracks] = energy_arr[selected_hit_idx];
        for (int d = 0; d < 3; d++) {
            event_info.momentum[processed_tracks][d] = momentum_arr[selected_hit_idx][d];
            event_info.track_origin[processed_tracks][d] = track_origin_arr[selected_hit_idx][d];
        }
        event_info.trigger_node[processed_tracks] = trigger_node_arr[selected_hit_idx];
        event_info.particle_id[processed_tracks] = particle_id_arr[selected_hit_idx];
        event_info.particle_type[processed_tracks] = particle_type_arr[selected_hit_idx];
        event_info.parent_particle_type[processed_tracks] = parent_particle_type_arr[selected_hit_idx];

        processed_tracks++;
    }
    event_info.num_tracks = processed_tracks;
    for (int d = 0; d < 3; d++) {
        event_info.interaction_point[d] = interaction_point_arr[d];
    }
    event_info.trigger = trigger;
    event_info.has_trigger_pair = has_trigger_pair;
}
