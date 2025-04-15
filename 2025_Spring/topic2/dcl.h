#ifndef DCL_H
#define DCL_H

#include <ap_fixed.h>
#include <hls_math.h>

// You may choose to use ap_fixed<16,5> if desired, but here we use float.
typedef float data_t;

// Define maximum sizes (adjust as necessary)
#define MAX_HITS       1024    // Maximum number of hits per event
#define MAX_EDGES      2048    // Maximum number of edges per event
#define MAX_TRACKS     256     // Maximum number of tracks (connected components)
#define MAX_TRACK_SIZE MAX_HITS  // Maximum number of hits in one track
#define NUM_LAYERS     5       // Number of detector layer groups

// Define the layer groups (matching Python's layers = [(0,), (1,), (2,), (3,4), (5,6)])
const int layer_start[NUM_LAYERS] = {0, 1, 2, 3, 5};
const int layer_end[NUM_LAYERS]   = {0, 1, 2, 4, 6};

// The EventInfo structure holds the processed event (track) attributes.
struct EventInfo {
    data_t n_pixels[MAX_TRACKS][NUM_LAYERS];       // Sum of pixel counts per track per layer group
    data_t energy[MAX_TRACKS];                     // Energy for each track
    data_t momentum[MAX_TRACKS][3];                // Momentum vector for each track
    data_t interaction_point[3];                   // Global interaction (collision) point
    bool trigger;                                  // Event trigger flag
    bool has_trigger_pair;                         // Indicates if the event has a paired trigger
    data_t track_origin[MAX_TRACKS][3];            // Estimated origin for each track
    int trigger_node[MAX_TRACKS];                  // Index of the trigger hit for each track
    int particle_id[MAX_TRACKS];                   // Representative particle ID for each track
    int particle_type[MAX_TRACKS];                 // Particle type classification for each track
    int parent_particle_type[MAX_TRACKS];          // Parent particle type for each track
    data_t track_hits[MAX_TRACKS][3 * NUM_LAYERS];   // Weighted hit positions (per layer group)
    int track_n_hits[MAX_TRACKS][NUM_LAYERS];        // Count of hits per layer group for each track
    int num_tracks;                                // Total number of tracks reconstructed
};

// Prototype for the HLS kernel function that reconstructs the tracks.
void compute_tracks_HLS(
    // Input arrays (all pre‑converted from NPZ to binary):
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
    // Output (reconstructed event information):
    EventInfo &event_info
);

#endif
