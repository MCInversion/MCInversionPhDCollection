import matplotlib.pyplot as plt
import numpy as np
import glob

# Define the file paths
file_paths = [
    "./BunnyLSW150FullWrap_Histograms.txt",
    "./BunnyLSW150Obstacle_Histograms.txt",
    "./BunnyRepulsive220_Histograms.txt",
    "./BunnyLSWDaniel_Histograms.txt"
]

# Define conversion factors
conversion_factors = {
    "Bunny Repulsive220": 1/10,
    "Bunny LSW200 Daniel": 1/1000,
    # Assuming no conversion needed for the reference datasets
    "Bunny LSW150 Full Wrap": 1, 
    "Bunny LSW150 Obstacle": 1
}

# Prepare data storage
hist_data = {
    "2 * inradius / circumradius": {},
    "Equilateral tri Jacobian condition number": {},
    "Vertex distance eval": {}
}

# Parse the files
for file_path in file_paths:
    with open(file_path, 'r') as file:
        lines = file.readlines()
        dataset_name = lines[0].strip().split(":")[1].strip()
        current_metric = ""
        for line in lines[1:]:
            if "Metric:" in line:
                current_metric = line.strip().split(":")[1].strip()
                hist_data[current_metric][dataset_name] = []
            elif line.startswith("[") or line.startswith("("):
                parts = line.split(":")
                range_part = parts[0].strip().strip("()[]").split("..")
                min_val, max_val = float(range_part[0]), float(range_part[1])
                count = int(parts[1].strip())
                if not np.isinf(min_val) and not np.isinf(max_val):
                    # Apply conversion factor if metric is "Vertex distance eval"
                    if current_metric == "Vertex distance eval":
                        conversion_factor = conversion_factors[dataset_name]
                        min_val *= conversion_factor
                        max_val *= conversion_factor
                    hist_data[current_metric][dataset_name].append((min_val, max_val, count))

# Create a figure for the plots
fig, axs = plt.subplots(len(hist_data), 1, figsize=(7, 5.5))

metric_index = 0
jacobian_x_cutoff = 8
for metric, datasets in hist_data.items():
    ax = axs[metric_index]
    for dataset_name, bins in datasets.items():
        # Adjust bin edges and counts for general case
        bin_edges = [bins[0][0]] + [b[0] for b in bins] + [bins[-1][1]]
        counts = [0] + [b[2] for b in bins] + [0]
        
        # Clip data at x=jacobian_x_cutoff for "Bunny LSW200 Daniel" in "Equilateral tri Jacobian condition number" metric
        if metric == "Equilateral tri Jacobian condition number" and dataset_name == "Bunny LSW200 Daniel":
            # Find the index where bin_edges exceed jacobian_x_cutoff and clip arrays
            clip_index = next((i for i, edge in enumerate(bin_edges) if edge > jacobian_x_cutoff), len(bin_edges))
            bin_edges = bin_edges[:clip_index+1]  # Include the boundary bin
            counts = counts[:clip_index+1]
            # Ensure the last visible bin edge is exactly at jacobian_x_cutoff
            bin_edges[-1] = jacobian_x_cutoff

        # Plot the histogram as a full curve
        ax.plot(bin_edges, counts, label=dataset_name, drawstyle='steps-pre')
        ax.fill_between(bin_edges, counts, step="pre", alpha=0.4)

    if metric == "Vertex distance eval":
        ax.set_ylabel('# of vertices')
    else:
        ax.set_ylabel('# of faces')

    ax.set_title(metric)
    ax.legend()
    ax.grid(True)
    metric_index += 1

plt.tight_layout()
plt.show()

