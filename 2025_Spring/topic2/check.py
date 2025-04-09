import numpy as np

data = np.load("./input_graphs/nontrigger/0/event000041188.npz")
print(data.files)  # ['edge_index', 'layer_id', 'n_pixels', 'hit_cartesian', ...]
