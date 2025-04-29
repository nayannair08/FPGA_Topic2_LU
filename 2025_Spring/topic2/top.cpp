#include "dcl.h"
// #define DEBUG

// Define a maximum number of neighbors per node.
#define MAX_NEIGHBORS 64

void filtering_edges(
    const int event,
    const int edge_index[NUM_EVENTS][2][MAX_EDGES],
    const data_t model_edge_probability[NUM_EVENTS][MAX_EDGES],
    const int num_edges[NUM_EVENTS],
    int edge_index_local[2][MAX_EDGES],
    data_t model_edge_probability_local[MAX_EDGES],
    int num_edges_local,
    int filtered_edges[2][MAX_EDGES],
    int filtered_count,
    int event1
) {
    event1 = event;
    //============================================================================================================
	//num_edges
	//edge_index
	//model_edge_probability
	//============================================================================================================
    // ----- Step 1: Filter Edges -----
    for (int i = 0; i < num_edges[event]; i++) {
        edge_index_local[0][i] = edge_index[event][0][i];
        edge_index_local[1][i] = edge_index[event][1][i];
        model_edge_probability_local[i] = model_edge_probability[event][i];
    }
    num_edges_local = num_edges[event];
    const data_t threshold = 0.5;
    for (int i = 0; i < num_edges_local; i++) {
        if (model_edge_probability_local[i] >= threshold) {
            filtered_edges[0][filtered_count] = edge_index_local[0][i];
            filtered_edges[1][filtered_count] = edge_index_local[1][i];
            filtered_count++;
        }
    }
}

void adjacency_list(
    const int num_nodes,
    const int filtered_edges[2][MAX_EDGES],
    const int filtered_count,
    int adj_list[MAX_HITS][MAX_NEIGHBORS],
    int adj_count[MAX_HITS],
    int max_node
) {
    //============================================================================================================
	//event
	//filtered edges
	//filtered count
	//
	//num_hits
	//============================================================================================================
    // ----- Step 2: Build the Adjacency List -----
    for (int i = 0; i < num_nodes; i++) {
        adj_count[i] = 0;
    }
    for (int i = 0; i < filtered_count; i++) {
        int u = filtered_edges[0][i];
        int v = filtered_edges[1][i];
        if (u > max_node) max_node = u;
        if (v > max_node) max_node = v;
        if (adj_count[u] < MAX_NEIGHBORS) {
            adj_list[u][adj_count[u]] = v;
            adj_count[u]++;
        }
        if (adj_count[v] < MAX_NEIGHBORS) {
            adj_list[v][adj_count[v]] = u;
            adj_count[v]++;
        }
    }
}

void iterative_dfs (
    const int num_nodes,
    const int max_node,
    const int adj_list[MAX_HITS][MAX_NEIGHBORS],
    const int adj_count[MAX_HITS],
    int tracks[MAX_TRACKS][MAX_TRACK_SIZE],
    int track_sizes[MAX_TRACKS],
    int track_count
) {
    //============================================================================================================
	//Adj list
	//Adj count
	//event
	//max_node
	//
	//============================================================================================================
    // ----- Step 3: Connected Components via Iterative DFS -----
    bool visited[MAX_HITS];
    for (int i = 0; i < num_nodes; i++) {
        visited[i] = false;
    }
    for (int i = 0; i < max_node+1; i++) {
        if (!visited[i]) {
            int current_track[MAX_TRACK_SIZE];
            int current_size = 0;
            int stack[MAX_HITS];
            int stack_ptr = 0;
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
            track_sizes[track_count] = current_size;
            track_count++;
        }
    }
}

void process_tracks (
    int event,
    const int track_count,
    const int tracks[MAX_TRACKS][MAX_TRACK_SIZE],
    const int track_sizes[MAX_TRACKS],
    const int layer_id[NUM_EVENTS][MAX_HITS],
    const int n_pixels_arr[NUM_EVENTS][MAX_HITS],
    const data_t hit_cartesian[NUM_EVENTS][MAX_HITS][3],
    const data_t particle_id_arr[NUM_EVENTS][MAX_HITS],
    const data_t energy_arr[NUM_EVENTS][MAX_HITS],
    const data_t momentum_arr[NUM_EVENTS][MAX_HITS][3],
    const data_t track_origin_arr[NUM_EVENTS][MAX_HITS][3],
    const int trigger_node_arr[NUM_EVENTS][MAX_HITS],
    const data_t particle_type_arr[NUM_EVENTS][MAX_HITS],
    const data_t parent_particle_type_arr[NUM_EVENTS][MAX_HITS],
    const data_t interaction_point_arr[NUM_EVENTS][3],
    const bool trigger[NUM_EVENTS],
    const bool has_trigger_pair[NUM_EVENTS],
    const bool intt_required[NUM_EVENTS],
    EventInfo event_info[NUM_EVENTS],
    int layer_id_local[MAX_HITS],
    int n_pixels_arr_local[MAX_HITS],
    data_t hit_cartesian_local[MAX_HITS][3],
    data_t particle_id_arr_local[MAX_HITS],
    data_t energy_arr_local[MAX_HITS],
    data_t momentum_arr_local[MAX_HITS][3],
    data_t track_origin_arr_local[MAX_HITS][3],
    int trigger_node_arr_local[MAX_HITS],
    data_t particle_type_arr_local[MAX_HITS],
    data_t parent_particle_type_arr_local[MAX_HITS],
    data_t interaction_point_arr_local[3],
    bool trigger_local,
    bool has_trigger_pair_local,
    bool intt_required_local
) {
    //============================================================================================================
	//event
	//track count
	//tracks
	//track sizes
	//
	//event_info
	//intt_required
	//layer_id
	//n_pixels_arr
	//hit_cartisian
	//particle_id_arr
	//energy_arr
	//momentum_arr
	//track_origin_arr
	//trigger_node_arr
	//particle_type_arr
	//parent_particle_type_arr
	//interaction_point_arr
	//trigger
	//has_trigger_pair
	//============================================================================================================
    for (int i = 0; i < MAX_HITS; i++) {
        layer_id_local[i] = layer_id[event][i];
        n_pixels_arr_local[i] = n_pixels_arr[event][i];
        for (int j = 0; j < 3; j++) {
            hit_cartesian_local[i][j] = hit_cartesian[event][i][j];
        }
        particle_id_arr_local[i] = particle_id_arr[event][i];
        energy_arr_local[i] = energy_arr[event][i];
        for (int j = 0; j < 3; j++) {
            momentum_arr_local[i][j] = momentum_arr[event][i][j];
            track_origin_arr_local[i][j] = track_origin_arr[event][i][j];
        }
        trigger_node_arr_local[i] = trigger_node_arr[event][i];
        particle_type_arr_local[i] = particle_type_arr[event][i];
        parent_particle_type_arr_local[i] = parent_particle_type_arr[event][i];
    }
    for (int i = 0; i < 3; i++) {
        interaction_point_arr_local[i] = interaction_point_arr[event][i];
    }
    trigger_local = trigger[event];
    has_trigger_pair_local = has_trigger_pair[event];
    intt_required_local = intt_required[event];
    // ----- Step 4: Process Each Track -----
    for (int i = 0; i < MAX_TRACKS; i++) {
        for (int j = 0; j < NUM_LAYERS; j++) {
            event_info[event].n_pixels[i][j] = 0;
            event_info[event].track_n_hits[i][j] = 0;
            for (int k = 0; k < 3; k++) {
                event_info[event].track_hits[i][3*j + k] = 0;
            }
        }
        event_info[event].energy[i] = 0;
        for (int k = 0; k < 3; k++) {
            event_info[event].momentum[i][k] = 0;
            event_info[event].track_origin[i][k] = 0;
        }
        event_info[event].trigger_node[i] = 0;
        event_info[event].particle_id[i] = 0;
        event_info[event].particle_type[i] = 0;
        event_info[event].parent_particle_type[i] = 0;
    }
    int processed_tracks = 0;
    for (int t = 0; t < track_count && processed_tracks < MAX_TRACKS; t++) {
        if (intt_required_local) {
            bool valid_track = false;
            for (int k = 0; k < track_sizes[t]; k++) {
                int hit_idx = tracks[t][k];
                if (layer_id_local[hit_idx] >= 3) {
                    valid_track = true;
                    break;
                }
            }
            if (!valid_track) continue;
        }
        for (int j = 0; j < NUM_LAYERS; j++) {
            data_t weighted_sum[3] = {0, 0, 0};
            int sum_pixels = 0;
            int count_hits = 0;
            for (int k = 0; k < track_sizes[t]; k++) {
                int hit_idx = tracks[t][k];
                int hit_layer = layer_id_local[hit_idx];
                if (hit_layer >= layer_start[j] && hit_layer <= layer_end[j]) {
                    int pix = n_pixels_arr_local[hit_idx];
                    sum_pixels += pix;
                    for (int d = 0; d < 3; d++) {
                        weighted_sum[d] += pix * hit_cartesian_local[hit_idx][d];
                    }
                    count_hits++;
                }
            }
            for (int d = 0; d < 3; d++) {
                if (sum_pixels > 0)
                    event_info[event].track_hits[processed_tracks][3*j + d] = weighted_sum[d] / sum_pixels;
                else
                    event_info[event].track_hits[processed_tracks][3*j + d] = 0;
            }
            event_info[event].n_pixels[processed_tracks][j] = sum_pixels;
            event_info[event].track_n_hits[processed_tracks][j] = count_hits;
        }
        data_t pids[MAX_TRACK_SIZE][3];
        for (int i = 0; i < MAX_TRACK_SIZE; i++) {
            pids[i][0] = -1;
            pids[i][1] = -1;
            pids[i][2] = -1;
        }
        for (int i = 0; i < track_sizes[t]; i++) {
            data_t p_id = particle_id_arr_local[tracks[t][i]];
            for (int j = 0; j < MAX_TRACK_SIZE; j++) {
                if (pids[j][0] == -1) {
                    pids[j][0] = p_id;
                    pids[j][1] = 1;
                    pids[j][2] = tracks[t][i];
                    break;
                } else if (pids[j][0] == p_id) {
                    pids[j][1]++;
                    break;
                }
            }
        }
        int max_count = 0;
        int mode_particle_id = -1;
        int track_particle_idx = -1;
        for (int i = 0; i < MAX_TRACK_SIZE; i++) {
            if (pids[i][1] > max_count) {
                max_count = pids[i][1];
                mode_particle_id = pids[i][0];
                track_particle_idx = pids[i][2];
            }
        }
        int rep_hit = tracks[t][track_particle_idx];
        event_info[event].energy[processed_tracks] = energy_arr_local[rep_hit];
        for (int d = 0; d < 3; d++) {
            event_info[event].momentum[processed_tracks][d] = momentum_arr_local[rep_hit][d];
            event_info[event].track_origin[processed_tracks][d] = track_origin_arr_local[rep_hit][d];
        }
        event_info[event].trigger_node[processed_tracks] = trigger_node_arr_local[rep_hit];
        event_info[event].particle_id[processed_tracks] = particle_id_arr_local[rep_hit];
        event_info[event].particle_type[processed_tracks] = particle_type_arr_local[rep_hit];
        event_info[event].parent_particle_type[processed_tracks] = parent_particle_type_arr_local[rep_hit];
        processed_tracks++;
    }
    event_info[event].num_tracks = processed_tracks;
    // ----- Step 5: Store Global Event Information -----
    for (int d = 0; d < 3; d++) {
        event_info[event].interaction_point[d] = interaction_point_arr_local[d];
    }
    event_info[event].trigger = trigger_local;
    event_info[event].has_trigger_pair = has_trigger_pair_local;
}




// ------------------------------
// Kernel Function: compute_tracks_HLS
// ------------------------------
// Reconstructs tracks from raw detector data using iterative DFS.
// "intt_required" indicates that only tracks with at least one hit with layer_id >= 3 are kept.
void compute_tracks_HLS(
    const int edge_index[NUM_EVENTS][2][MAX_EDGES],
    const data_t model_edge_probability[NUM_EVENTS][MAX_EDGES],
    const int num_edges[NUM_EVENTS],
    const int layer_id[NUM_EVENTS][MAX_HITS],
    const int n_pixels_arr[NUM_EVENTS][MAX_HITS],
    const data_t hit_cartesian[NUM_EVENTS][MAX_HITS][3],
    const data_t particle_id_arr[NUM_EVENTS][MAX_HITS],
    const data_t energy_arr[NUM_EVENTS][MAX_HITS],
    const data_t momentum_arr[NUM_EVENTS][MAX_HITS][3],
    const data_t track_origin_arr[NUM_EVENTS][MAX_HITS][3],
    const int trigger_node_arr[NUM_EVENTS][MAX_HITS],
    const data_t particle_type_arr[NUM_EVENTS][MAX_HITS],
    const data_t parent_particle_type_arr[NUM_EVENTS][MAX_HITS],
    const int num_hits[NUM_EVENTS],
    const data_t interaction_point_arr[NUM_EVENTS][3],
    bool trigger[NUM_EVENTS],
    bool has_trigger_pair[NUM_EVENTS],
    bool intt_required[NUM_EVENTS],
    EventInfo event_info[NUM_EVENTS]
) {

	//ToDo: Interface Pragmas
	
#pragma HLS interface m_axi port=edge_index offset=slave bundle=mem1
#pragma HLS interface m_axi port=model_edge_probability offset=slave bundle=mem1
#pragma HLS interface m_axi port=num_edges offset=slave bundle=mem1

#pragma HLS interface m_axi port=num_hits offset=slave bundle=mem2

#pragma HLS interface m_axi port=layer_id offset=slave bundle=mem3
#pragma HLS interface m_axi port=n_pixels_arr offset=slave bundle=mem3
#pragma HLS interface m_axi port=hit_cartesian offset=slave bundle=mem3
#pragma HLS interface m_axi port=particle_id_arr offset=slave bundle=mem3
#pragma HLS interface m_axi port=energy_arr offset=slave bundle=mem3
#pragma HLS interface m_axi port=momentum_arr offset=slave bundle=mem3
#pragma HLS interface m_axi port=track_origin_arr offset=slave bundle=mem3
#pragma HLS interface m_axi port=trigger_node_arr offset=slave bundle=mem3
#pragma HLS interface m_axi port=particle_type_arr offset=slave bundle=mem3
#pragma HLS interface m_axi port=parent_particle_type_arr offset=slave bundle=mem3
#pragma HLS interface m_axi port=interaction_point_arr offset=slave bundle=mem3
#pragma HLS interface m_axi port=trigger offset=slave bundle=mem3
#pragma HLS interface m_axi port=has_trigger_pair offset=slave bundle=mem3
#pragma HLS interface m_axi port=intt_required offset=slave bundle=mem3
#pragma HLS interface m_axi port=event_info offset=slave bundle=mem3

#pragma HLS interface s_axilite port=return

    //Local Data Structures
    int edge_index_local[2][MAX_EDGES];
    data_t model_edge_probability_local[MAX_EDGES];
    int num_edges_local = 0;
    int layer_id_local[MAX_HITS];
    int n_pixels_arr_local[MAX_HITS];
    data_t hit_cartesian_local[MAX_HITS][3];
    data_t particle_id_arr_local[MAX_HITS];
    data_t energy_arr_local[MAX_HITS];
    data_t momentum_arr_local[MAX_HITS][3];
    data_t track_origin_arr_local[MAX_HITS][3];
    int trigger_node_arr_local[MAX_HITS];
    data_t particle_type_arr_local[MAX_HITS];
    data_t parent_particle_type_arr_local[MAX_HITS];
    int num_hits_local = 0;
    data_t interaction_point_arr_local[3];
    bool trigger_local = false;
    bool has_trigger_pair_local = false;
    bool intt_required_local = false;


    for (int event = 0; event < NUM_EVENTS; event++) {
		
        int filtered_edges[2][MAX_EDGES];
        int filtered_count = 0;
        int adj_list[MAX_HITS][MAX_NEIGHBORS];
        int adj_count[MAX_HITS];
        num_hits_local = num_hits[event];
        int num_nodes = num_hits_local; // Each hit is a node.
        int max_node = 0;
        int tracks[MAX_TRACKS][MAX_TRACK_SIZE];
        int track_sizes[MAX_TRACKS];
        int track_count = 0;
        int event1 = 0;

        // ----- Step 1: Filter Edges -----
        filtering_edges(
            event,
            edge_index,
            model_edge_probability,
            num_edges,
            edge_index_local,
            model_edge_probability_local,
            num_edges_local,
            filtered_edges,
            filtered_count,
            event1
        );
        // ----- Step 2: Build the Adjacency List -----
        adjacency_list(
            num_nodes,
            filtered_edges,
            filtered_count,
            adj_list,
            adj_count,
            max_node
        );
        // ----- Step 3: Connected Components via Iterative DFS -----
        iterative_dfs(
            num_nodes,
            max_node,
            adj_list,
            adj_count,
            tracks,
            track_sizes,
            track_count
        );
        // ----- Step 4: Process Each Track -----
        process_tracks(
            event1,
            track_count,
            tracks,
            track_sizes,
            layer_id,
            n_pixels_arr,
            hit_cartesian,
            particle_id_arr,
            energy_arr,
            momentum_arr,
            track_origin_arr,
            trigger_node_arr,
            particle_type_arr,
            parent_particle_type_arr,
            interaction_point_arr,
            trigger,
            has_trigger_pair,
            intt_required,
            event_info,
            layer_id_local,
            n_pixels_arr_local,
            hit_cartesian_local,
            particle_id_arr_local,
            energy_arr_local,
            momentum_arr_local,
            track_origin_arr_local,
            trigger_node_arr_local,
            particle_type_arr_local,
            parent_particle_type_arr_local,
            interaction_point_arr_local,
            trigger_local,
            has_trigger_pair_local,
            intt_required_local
        );
    }

}