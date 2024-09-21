import os
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
import imageio.v2 as imageio  # Use imageio.v2 to avoid the deprecation warning
from collections import defaultdict
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg

# Define the procedure name and the path to your directory containing the files
#procedure_name = "ShrinkWrappingAnIncompleteCirclePointCloud_NoInnerCircleNoRemeshing"
#procedure_name = "ShrinkWrappingAnIncompleteCirclePointCloud_NoRemeshing"
#procedure_name = "ShrinkingAndExpandingCircle_NoRemeshing"
#procedure_name = "ShrinkingAndExpandingCircle_SlowerInnerEtaNoRemeshing"
#procedure_name = "ShrinkWrappingAnIncompleteCirclePointCloud_WithRemeshing"
#procedure_name = "ShrinkWrappingAnIncompleteDeformedCirclePointCloud_WithRemeshing"
#procedure_name = "InnerOuterLSWOfImportedMeshPtCloudSlice_WithRemeshing"
#procedure_name = "armadillonewEvol_Pts2D3"
#procedure_name = "bunnynewEvol_Pts2D3"
#procedure_name = "maxPlancknewEvol_Pts2D3"
procedure_name = "nefertitinewEvol_Pts2D3"

#directory = "../output/core_tests/"  # Adjust this path accordingly
directory = "../output"  # Adjust this path accordingly

# Helper function to read the grid dimensions from the *.gdim2d file
def read_gdim2d_file(filepath):
    bbox_min = bbox_max = (0, 0)
    nx = ny = 0
    cell_size = 0.0

    with open(filepath, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if "Grid Dimensions" in line:
                parts = line.split(':')[1].strip().split(',')
                nx = int(parts[0].strip())
                ny = int(parts[1].strip())
            elif "Bounding Box" in line:
                parts = line.split(':')[1].strip().split('->')
                bbox_min = tuple(map(float, parts[0].strip('() ').split(',')))
                bbox_max = tuple(map(float, parts[1].strip('() ').split(',')))
            elif "Cell Size" in line:
                cell_size = float(line.split(':')[1].strip())

    return bbox_min, bbox_max, nx, ny, cell_size

# Load the background image and grid dimensions if both files exist
background_image_path = os.path.join(directory, f"{procedure_name}_TargetDF.png")
legend_image_path = os.path.join(directory, f"{procedure_name}_TargetDF_Scale.png")
gdim2d_path = os.path.join(directory, f"{procedure_name}_TargetDF.gdim2d")

if os.path.exists(background_image_path) and os.path.exists(gdim2d_path) and os.path.exists(legend_image_path):
    bbox_min, bbox_max, nx, ny, cell_size = read_gdim2d_file(gdim2d_path)
    background_img = imageio.imread(background_image_path)
    extent = [bbox_min[0], bbox_max[0], bbox_max[1], bbox_min[1]]

    legend_img = imageio.imread(legend_image_path)
    bbox_size = np.array(bbox_max) - np.array(bbox_min)
    clip_factor = 0.05
    legend_width = legend_img.shape[1]  # Get the image width in pixels
    clipped_width = int(legend_width * clip_factor)
    legend_img = legend_img[:, clipped_width:] # clip legend
    extent_legend = [bbox_max[0] - 0.175 * bbox_size[1],
                     bbox_max[0],
                     bbox_min[1] + 0.15 * bbox_size[1],
                     bbox_max[1] - 0.15 * bbox_size[1]]
else:
    background_img = None

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

# Gather all the PLY files for inner and outer curves in sorted order
inner_ply_files = defaultdict(dict)  # Use a dictionary of dictionaries
time_ids = []  # To store the time_id values

for f in sorted(os.listdir(directory)):
    if f.startswith(f"{procedure_name}_Inner") and f.endswith(".ply"):
        # Extract the curve ID and time ID from the filename
        parts = f.split('_')
        curve_id = parts[2].replace('Inner', '')  # Extract curve ID (e.g., "0" from "Inner0")
        time_id = int(parts[-1].replace('Evol', '').replace('.ply', ''))  # Extract time ID (e.g., "1" from "Evol1")
        inner_ply_files[curve_id][time_id] = os.path.join(directory, f)
        if time_id not in time_ids:
            time_ids.append(time_id)  # Collect time_id if not already added

outer_ply_files = defaultdict(str)  # Use a dictionary for outer files indexed by time_id

for f in sorted(os.listdir(directory)):
    if f.startswith(f"{procedure_name}_Outer_Evol") and f.endswith(".ply"):
        # Extract the time ID from the filename
        time_id = int(f.split('_')[-1].replace('Evol', '').replace('.ply', ''))  # Extract time ID
        outer_ply_files[time_id] = os.path.join(directory, f)
        if time_id not in time_ids:
            time_ids.append(time_id)  # Collect time_id if not already added

# Sort inner files by time_id within each curve_id
for curve_id in inner_ply_files:
    inner_ply_files[curve_id] = dict(sorted(inner_ply_files[curve_id].items()))

# Sort outer files by time_id
outer_ply_files = dict(sorted(outer_ply_files.items()))

# Sort the time_ids
time_ids = sorted(time_ids)

# Read all vertex positions from the PLY files for inner and outer curves
inner_curves = {curve_id: [read_ply(ply_file) for ply_file in files.values()]
                for curve_id, files in inner_ply_files.items()}
outer_curves = [read_ply(ply_file) for ply_file in outer_ply_files.values()]

# Debugging: Print the number of vertices read from each file
for inner_id, inner_curve_group in inner_curves.items():
    for i, inner_curve in enumerate(inner_curve_group):
        print(f"Inner file [{inner_id}][{i}]: {inner_ply_files[inner_id][i]}: {len(inner_curve)} vertices")

for i, outer_curve in enumerate(outer_curves):
    print(f"Outer file [{i}]: {outer_ply_files[i]}: {len(outer_curve)} vertices")

# Set up the figure with two subplots: one for the main plot, one for the legend
fig = plt.figure(figsize=(8, 6))
gs = gridspec.GridSpec(1, 2, width_ratios=[8, 1])  # ratio for main plot and legend
plt.subplots_adjust(wspace=-0.3)

# Create the main plot in the first grid section
ax_main = plt.subplot(gs[0])

# Create a separate axis for the legend in the second grid section
ax_legend = plt.subplot(gs[1])

# Disable the axes for the legend (since it's an image)
ax_legend.axis('off')

# Display the background image if it was successfully loaded
if background_img is not None:
    ax_main.imshow(background_img, extent=extent)
    ax_main.set_xlim(extent[:2])
    ax_main.set_ylim(extent[2:])
    ax_main.invert_yaxis()  # Reverse the direction of the y-axis
else:
    print("missing background_img. Resizing automatically according to curves.")
    # Automatically adjust the plot limits based on the curves if no background image
    all_points = np.concatenate([curve for curve in outer_curves] + 
                                [curve for curves in inner_curves.values() for curve in curves])
    ax_main.set_xlim(np.min(all_points[:, 0]), np.max(all_points[:, 0]))
    ax_main.set_ylim(np.min(all_points[:, 1]), np.max(all_points[:, 1]))
    ax_main.invert_yaxis()  # Reverse the direction of the y-axis

# Show legend
if legend_img is not None:
    ax_legend.imshow(legend_img, extent=extent_legend)

    # Load and overlay the symbol above the legend within the same axis
    legend_symbol = "./dGamma_200dpi.png"  # Path to your symbol file
    symbol_img = mpimg.imread(legend_symbol)
    
    # Overlay the symbol in the upper portion of the same axis
    symbol_extent = [extent_legend[0] + 0.04 * bbox_size[0], 
                     extent_legend[1] - 0.07 * bbox_size[0], 
                     extent_legend[3] + 0.02 * bbox_size[1],
                     extent_legend[3] + 0.1 * bbox_size[1]
    ]  # Adjust these values to control position
    ax_legend.imshow(legend_img, extent=extent_legend, aspect='equal')
    ax_legend.imshow(symbol_img, extent=symbol_extent, aspect='equal')
    ax_legend.set_xlim(extent_legend[:2])
    ax_legend.set_ylim([extent_legend[2], extent_legend[3] + 0.1 * bbox_size[1]])

ax_main.set_aspect('equal')

# Initialize the lines for inner (red) and outer (black) curves
outer_line, = ax_main.plot([], [], 'k-', linewidth=2)  # Black outer curve
inner_line_handles = [ax_main.plot([], [], '-', linewidth=2)[0] for _ in inner_curves]  # Different lines for each inner curve

# Initialize the line data
def init():
    outer_line.set_data([], [])
    for inner_lines in inner_line_handles:
        inner_lines.set_data([], [])
    return [outer_line] + inner_line_handles

# Update function for animation
def update(frame):
    # Update outer curve
    if outer_curves[frame].size > 0:
        x_outer, y_outer = outer_curves[frame].T
        outer_line.set_data(x_outer, y_outer)
    
    # Update inner curves
    for inner_lines, inner_curve_group in zip(inner_line_handles, inner_curves.values()):
        if inner_curve_group[frame].size > 0:
            x_inner, y_inner = inner_curve_group[frame].T
            inner_lines.set_data(x_inner, y_inner)
    
    # Update the title with the current time step using time_id
    ax_main.set_title(f"Time Step: {time_ids[frame]}", fontsize=14)

    return [outer_line] + inner_line_handles

# Create the animation
ani = FuncAnimation(fig, update, frames=len(outer_curves), init_func=init, blit=False, repeat=True)

# Save the animation as a GIF
output_gif_path = os.path.join(directory, f"{procedure_name}_animation.gif")
ani.save(output_gif_path, writer='pillow', fps=30)

# Show the animation
plt.show()
