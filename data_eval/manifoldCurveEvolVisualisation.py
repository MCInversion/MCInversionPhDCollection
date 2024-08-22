import os
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
import imageio.v2 as imageio  # Use imageio.v2 to avoid the deprecation warning
from collections import defaultdict

# Define the procedure name and the path to your directory containing the files
#procedure_name = "ShrinkWrappingAnIncompleteCirclePointCloud_NoInnerCircleNoRemeshing"
procedure_name = "ShrinkWrappingAnIncompleteCirclePointCloud_NoRemeshing"
directory = "../output/core_tests/"  # Adjust this path accordingly

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
gdim2d_path = os.path.join(directory, f"{procedure_name}_TargetDF.gdim2d")

if os.path.exists(background_image_path) and os.path.exists(gdim2d_path):
    bbox_min, bbox_max, nx, ny, cell_size = read_gdim2d_file(gdim2d_path)
    background_img = imageio.imread(background_image_path)
    extent = [bbox_min[0], bbox_max[0], bbox_min[1], bbox_max[1]]
else:
    background_img = None

# Helper function to read PLY files and extract vertices
def read_ply(file_path):
    vertices = []
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()

        vertex_section = False
        for line in lines:
            if line.startswith('end_header'):
                vertex_section = True
                continue
            if vertex_section and len(line.strip()) > 0:
                coords = list(map(float, line.strip().split()))
                if len(coords) == 3:
                    vertices.append((coords[0], coords[1]))

        # Ensure the curve is closed by connecting the last vertex to the first
        if len(vertices) > 1 and vertices[-1] != vertices[0]:
            vertices.append(vertices[0])

    except Exception as e:
        print(f"Error reading {file_path}: {e}")
    return np.array(vertices)

# Gather all the PLY files for inner and outer curves in sorted order
inner_ply_files = defaultdict(dict)  # Use a dictionary of dictionaries

for f in sorted(os.listdir(directory)):
    if f.startswith(f"{procedure_name}_Inner") and f.endswith(".ply"):
        # Extract the curve ID and time ID from the filename
        parts = f.split('_')
        curve_id = parts[2].replace('Inner', '')  # Extract curve ID (e.g., "0" from "Inner0")
        time_id = int(parts[-1].replace('Evol', '').replace('.ply', ''))  # Extract time ID (e.g., "1" from "Evol1")
        inner_ply_files[curve_id][time_id] = os.path.join(directory, f)

outer_ply_files = defaultdict(str)  # Use a dictionary for outer files indexed by time_id

for f in sorted(os.listdir(directory)):
    if f.startswith(f"{procedure_name}_Outer_Evol") and f.endswith(".ply"):
        # Extract the time ID from the filename
        time_id = int(f.split('_')[-1].replace('Evol', '').replace('.ply', ''))  # Extract time ID
        outer_ply_files[time_id] = os.path.join(directory, f)

# Sort inner files by time_id within each curve_id
for curve_id in inner_ply_files:
    inner_ply_files[curve_id] = dict(sorted(inner_ply_files[curve_id].items()))

# Sort outer files by time_id
outer_ply_files = dict(sorted(outer_ply_files.items()))

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

# Set up the plot
fig, ax = plt.subplots()

# Display the background image if it was successfully loaded
if background_img is not None:
    ax.imshow(background_img, extent=extent)
    ax.set_xlim(extent[:2])
    ax.set_ylim(extent[2:])
else:
    print("missing background_img. Resizing automatically according to curves.")
    # Automatically adjust the plot limits based on the curves if no background image
    all_points = np.concatenate([curve for curve in outer_curves] + 
                                [curve for curves in inner_curves.values() for curve in curves])
    ax.set_xlim(np.min(all_points[:, 0]), np.max(all_points[:, 0]))
    ax.set_ylim(np.min(all_points[:, 1]), np.max(all_points[:, 1]))

ax.set_aspect('equal')

# Initialize the lines for inner (red) and outer (black) curves
outer_line, = ax.plot([], [], 'k-', linewidth=2)  # Black outer curve
inner_line_handles = [ax.plot([], [], '-', linewidth=2)[0] for _ in inner_curves]  # Different lines for each inner curve

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

    return [outer_line] + inner_line_handles

# Create the animation
ani = FuncAnimation(fig, update, frames=len(outer_curves), init_func=init, blit=True, repeat=True)

# Save the animation as a GIF
output_gif_path = os.path.join(directory, f"{procedure_name}_animation.gif")
ani.save(output_gif_path, writer='pillow', fps=2)

# Show the animation
plt.show()
