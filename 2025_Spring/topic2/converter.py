import numpy as np
import os

def write_array_to_bin(array, filename):
    array = np.asarray(array)
    
    # Convert booleans to uint8 for hardware compatibility
    if array.dtype == np.bool_:
        array = array.astype(np.uint8)

    # Flatten and write as float32 or int32 depending on dtype
    if np.issubdtype(array.dtype, np.integer):
        array = array.astype(np.int32)
    else:
        array = array.astype(np.float32)

    with open(filename, 'wb') as f:
        f.write(array.tobytes())

def convert_npz_to_bin(npz_path, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    data = np.load(npz_path, allow_pickle=False)

    for key in data.files:
        filename = os.path.join(out_dir, f"{key}.bin")
        write_array_to_bin(data[key], filename)
        print(f"Wrote {key} to {filename}")

# Example usage
if __name__ == '__main__':
    input_npz = "./input_graphs/nontrigger/0/event000041188.npz"  # replace with your input file
    output_dir = "event000041188"  # replace with desired output directory
    convert_npz_to_bin(input_npz, output_dir)
