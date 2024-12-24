import os
import matplotlib.pyplot as plt
import numpy as np
import imageio.v2 as imageio
from matplotlib.patches import Rectangle

draw_grid_lines = True  # Global flag to toggle grid lines
grid_color = "gray"  # Color for grid lines

def read_gdim2d_file(filepath):
    """Read grid dimensions from the .gdim2d file."""
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

        #print(edges)

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

# Read PLY file and extract only points
def read_ply_pts(file_path):
    """Read PLY file and extract only points."""
    vertices = []
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()

        num_vertices = 0

        # Parse header to get number of vertices
        for i, line in enumerate(lines):
            if line.startswith('element vertex'):
                num_vertices = int(line.split()[2])
            elif line.startswith('end_header'):
                header_end = i
                break

        # Read vertices
        for line in lines[header_end + 1:header_end + 1 + num_vertices]:
            coords = list(map(float, line.strip().split()))
            vertices.append((coords[0], coords[1]))

    except Exception as e:
        print(f"Error reading {file_path}: {e}")
    return np.array(vertices)

# Function to draw white grid lines
def draw_grid_lines(ax, bbox_min, bbox_max, nx, ny, cell_size):
    # Draw vertical lines
    for i in range(nx + 1):
        x = bbox_min[0] + i * cell_size
        ax.plot([x, x], [bbox_min[1], bbox_max[1]], color=grid_color, linewidth=0.5, alpha=0.5)

    # Draw horizontal lines
    for i in range(ny + 1):
        y = bbox_min[1] + i * cell_size
        ax.plot([bbox_min[0], bbox_max[0]], [y, y], color=grid_color, linewidth=0.5, alpha=0.5)

# Main visualization function
def visualize_field_and_curve(directory, base_filename):
    """Visualize the distance field and its corresponding curve."""
    gdim2d_path = os.path.join(directory, f"{base_filename}.gdim2d")
    distance_field_path = os.path.join(directory, f"{base_filename}.png")
    curve_path = os.path.join(directory, f"{base_filename}.ply")

    if not os.path.exists(gdim2d_path) or not os.path.exists(distance_field_path) or not os.path.exists(curve_path):
        print("Required files are missing.")
        return

    # Read grid dimensions and bounding box
    bbox_min, bbox_max, nx, ny, cell_size = read_gdim2d_file(gdim2d_path)
    # Adjust the extent using the real-world coordinates, scaled by pixel dimensions
    # The extent is calculated as bbox_min to (bbox_min + (number of cells * cell size))
    extent = [
        bbox_min[0],  # X-min
        bbox_min[0] + nx * cell_size,  # X-max
        bbox_min[1],  # Y-min
        bbox_min[1] + ny * cell_size  # Y-max
    ]

    # Load distance field image
    distance_field_img = imageio.imread(distance_field_path)
    distance_field_img = np.flipud(distance_field_img)

    # Determine if the PLY file contains edges or just points
    with open(curve_path, 'r') as f:
        contains_edges = any("element edge" in line for line in f)

    # Load curve or point data
    if contains_edges:
        curve_data = read_ply(curve_path)
        is_polyline = True
    else:
        curve_data = read_ply_pts(curve_path)
        is_polyline = False

    # Plot the field and curve or points
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.imshow(distance_field_img, extent=extent, origin='upper')

    if is_polyline:
        ax.plot(curve_data[:, 0], curve_data[:, 1], color='black', linewidth=2, label='Curve')
    else:
        ax.scatter(curve_data[:, 0], curve_data[:, 1], color='black', s=10, label='Point Cloud')

    ax.set_xlim([bbox_min[0], bbox_max[0]])
    ax.set_ylim([bbox_min[1], bbox_max[1]])

    # Draw grid lines if enabled
    if draw_grid_lines:
        draw_grid_lines(ax, bbox_min, bbox_max, nx, ny, cell_size)

    # Add legend and labels
    ax.legend()
    ax.set_title(f"Distance Field and Curve: {base_filename}")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")

    plt.tight_layout()

    # Save the plot as a PNG
    output_png_path = os.path.join(directory, f"{base_filename}_SDF.png")
    print("Saving ", output_png_path)
    plt.savefig(output_png_path, dpi=300)

    plt.show()

# Example usage
directory = "../output/sdf_tests"

#base_filename = "SimpleClosedBaseCurve"
#base_filename = "SimpleClosedManifoldCurve"
#base_filename = "SimpleOpenBaseCurve"
#base_filename = "SimpleOpenManifoldCurve"
#base_filename = "SimplePointCloud"
#base_filename = "SimpleClosedBaseCurve"
#base_filename = "SimpleClosedBaseCurveUnsigned"
#base_filename = "MoreComplexPolygonalManifoldCurve"
#base_filename = "MoreComplexPointCloud"
#base_filename = "SimpleClosedBaseCurveUnsignedInFirstQuadrant"
#base_filename = "SimpleClosedManifoldCurveUnsignedInFirstQuadrant"
base_filename = "ManifoldCircleCurve"
#base_filename = "CirclePointCloud"

visualize_field_and_curve(directory, base_filename)
