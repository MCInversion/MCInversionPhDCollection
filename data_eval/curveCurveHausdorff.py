import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import directed_hausdorff
from collections import defaultdict

# List of procedure names
procedure_names = [
    "curve0_Repulsionless",
    "curve1_Repulsionless",
    "equilibriumPair0",
    "equilibriumPair1",
    "equilibriumPair2",
    "equilibriumPair3"
]

# Map of procedure names to their proper names for the legend
procedure_name_map = {
    "curve0_Repulsionless": "(NEQ)",
    "curve1_Repulsionless": "(EQ1)",
    "equilibriumPair0": "(EQ2)",
    "equilibriumPair1": "(EQ4)",
    "equilibriumPair2": "(EQ3)",
    "equilibriumPair3": "(EQ5)"
}

# Gather all the PLY files for inner and outer curves in sorted order
directory = "../output"  # Adjust this path accordingly

# Helper function to read PLY files and extract vertices
def read_ply(file_path):
    vertices = []
    edges = []
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()

        num_vertices = 0
        num_edges = 0

        # Parse header to get number of vertices and edges
        for i, line in enumerate(lines):
            if line.startswith('element vertex'):
                num_vertices = int(line.split()[2])
            elif line.startswith('element edge'):
                num_edges = int(line.split()[2])
            elif line.startswith('end_header'):
                header_end = i
                break

        # Read vertices
        for line in lines[header_end + 1:header_end + 1 + num_vertices]:
            coords = list(map(float, line.strip().split()))
            vertices.append((coords[0], coords[1]))

        # Read edges
        for line in lines[header_end + 1 + num_vertices:header_end + 1 + num_vertices + num_edges]:
            edge = list(map(int, line.strip().split()))
            edges.append(edge)

        # Initialize polyline and visited sets
        polyline = []
        visited_vertices = set()
        visited_edges = set()

        # Start with the first edge
        current_vertex = edges[0][0]
        polyline.append(vertices[current_vertex])
        visited_vertices.add(current_vertex)

        # Build the polyline by following the edges in sequence
        while len(visited_edges) < len(edges):
            for i, edge in enumerate(edges):
                if i in visited_edges:
                    continue

                v0, v1 = edge[0], edge[1]

                # Follow the edge that connects to the current vertex
                if v0 == current_vertex and v1 not in visited_vertices:
                    polyline.append(vertices[v1])
                    visited_vertices.add(v1)
                    visited_edges.add(i)
                    current_vertex = v1
                    break
                elif v1 == current_vertex and v0 not in visited_vertices:
                    polyline.append(vertices[v0])
                    visited_vertices.add(v0)
                    visited_edges.add(i)
                    current_vertex = v0
                    break
            else:
                # If no edges were added, we've completed a loop or all vertices are visited
                break

        # Ensure the polyline is closed by adding the starting point at the end
        if polyline[-1] != polyline[0]:
            polyline.append(polyline[0])

    except Exception as e:
        print(f"Error reading {file_path}: {e}")
    return np.array(polyline)

# Function to compute the Hausdorff distance between two curves
def hausdorff_distance(curve1, curve2):
    def pointwise_distance(a, b):
        return np.linalg.norm(a - b, axis=1)

    dist_1_to_2 = np.max([np.min(pointwise_distance(pt, curve2)) for pt in curve1])
    dist_2_to_1 = np.max([np.min(pointwise_distance(pt, curve1)) for pt in curve2])
    return max(dist_1_to_2, dist_2_to_1)

# Function to compute Hausdorff distances for all time steps and plot the results
def compute_and_plot_hausdorff_distances(procedure_names):
    plt.figure(figsize=(5, 3))

    for procedure_name in procedure_names:
        print(f"Processing {procedure_name}...")

        # Initialize dictionaries for inner and outer curves
        inner_ply_files = defaultdict(dict)
        outer_ply_files = defaultdict(str)
        time_ids = []

        # Gather PLY files
        for f in sorted(os.listdir(directory)):
            if f.startswith(f"{procedure_name}_Inner") and f.endswith(".ply"):
                parts = f.split('_')
                curve_id_part = [part for part in parts if "Inner" in part][0]
                curve_id = curve_id_part.replace('Inner', '')  # Extract curve ID
                time_id = int(parts[-1].replace('Evol', '').replace('.ply', ''))  # Extract time ID
                inner_ply_files[curve_id][time_id] = os.path.join(directory, f)
                if time_id not in time_ids:
                    time_ids.append(time_id)

            if f.startswith(f"{procedure_name}_Outer_Evol") and f.endswith(".ply"):
                time_id = int(f.split('_')[-1].replace('Evol', '').replace('.ply', ''))  # Extract time ID
                outer_ply_files[time_id] = os.path.join(directory, f)
                if time_id not in time_ids:
                    time_ids.append(time_id)

        time_ids = sorted(time_ids)
        hausdorff_distances = []

        # Compute Hausdorff distances for each time step
        for time_id in time_ids:
            if time_id in outer_ply_files:
                outer_curve = read_ply(outer_ply_files[time_id])
            else:
                continue  # Skip if outer curve is missing

            inner_curve = None
            for curve_id, inner_files in inner_ply_files.items():
                if time_id in inner_files:
                    inner_curve = read_ply(inner_files[time_id])
                    break

            if inner_curve is not None and outer_curve is not None:
                h_dist = hausdorff_distance(outer_curve, inner_curve)
                hausdorff_distances.append(h_dist)

        # Normalize Hausdorff distances
        max_distance = max(hausdorff_distances) if hausdorff_distances else 1
        normalized_distances = [(d / max_distance) * 100 for d in hausdorff_distances]

        # Plot normalized Hausdorff distances for this procedure
        proper_name = procedure_name_map.get(procedure_name, procedure_name)
        plt.plot(time_ids[:len(normalized_distances)], normalized_distances, label=f'{proper_name}', linewidth=2)

    # Customize the plot
    #plt.title('Normalized Hausdorff Distance Over Time for Different Procedures')
    plt.xlabel('Time Step')
    plt.ylabel("$d_{{H}}$(F, G)")
    plt.legend(loc='upper center')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

# Call the function to compute and plot Hausdorff distances
compute_and_plot_hausdorff_distances(procedure_names)