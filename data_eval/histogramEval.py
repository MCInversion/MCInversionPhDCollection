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
fig, axs = plt.subplots(len(hist_data), 1, figsize=(8, 6))

metric_index = 0
for metric, datasets in hist_data.items():
    ax = axs[metric_index]
    for dataset_name, bins in datasets.items():
        bin_edges = [b[0] for b in bins] + [bins[-1][1]]
        counts = [b[2] for b in bins]
        # Adjust the label to include conversion factor if applicable
        label = f"{dataset_name}"
        ax.plot(bin_edges[:-1], counts, label=label, drawstyle='steps-post')
    
    # Here, we use 'metric' instead of 'current_metric'
    if metric == "Vertex distance eval":  # Make sure this matches exactly with your metric name
        ax.set_ylabel('Vertices')
    else:
        ax.set_ylabel('Faces')
    
    ax.set_title(metric)
    ax.legend()
    ax.grid(True)
    metric_index += 1

plt.tight_layout()
plt.show()

