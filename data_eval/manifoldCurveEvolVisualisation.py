import os
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
import imageio.v2 as imageio  # Use imageio.v2 to avoid the deprecation warning
from collections import defaultdict
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg
import matplotlib.patches as patches

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: Don't put keywords "Inner" or "Outer" in the procedure name because it's used to identify inner/outer curves
# ====================================================================================================

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
#procedure_name = "nefertitinewEvol_Pts2D3"

#procedure_name = "armadillonewEvol_Pts2D3_Repulsionless"
#procedure_name = "bunnynewEvol_Pts2D3_Repulsionless"
#procedure_name = "maxPlancknewEvol_Pts2D3_Repulsionless"

#procedure_name = "nefertitinewEvol_Pts2D3_Repulsionless"
#procedure_name = "nefertitinewEvol_Pts2D3_RepulsionlessTwo"

#procedure_name = "armadillonewEvol_Pts2D3_WithRepulsion"
#procedure_name = "bunnynewEvol_Pts2D3_WithRepulsion"
#procedure_name = "maxPlancknewEvol_Pts2D3_WithRepulsion"
#procedure_name = "nefertitinewEvol_Pts2D3_WithRepulsion"

#procedure_name = "curve0_Repulsionless"
#procedure_name = "curve1_Repulsionless"
#procedure_name = "curve2_Repulsionless"
#procedure_name = "curve3_Repulsionless"

#procedure_name = "concentricCircles0_Repulsionless"
#procedure_name = "concentricCircles1_Repulsionless"
#procedure_name = "concentricCircles2_Repulsionless"
procedure_name = "concentricCircles3_Repulsionless"

#procedure_name = "singleInnerCircleTest"
#procedure_name = "innerCircleOuterCirclePtsTest"

#procedure_name = "innerCircle_circle_PtsTest"
#procedure_name = "innerCircle_incompleteCircle_PtsTest"
#procedure_name = "innerCircle_sineDeformedCircle_PtsTest"
#procedure_name = "innerCircle_sineDeformedIncompleteCircle_PtsTest"
#procedure_name = "innerCircle_chamferedRectangle_PtsTest"
#procedure_name = "innerCircle_incompleteChamferedRectangle_PtsTest"
#procedure_name = "innerCircle_chamferedTriangle_PtsTest"
#procedure_name = "innerCircle_incompleteChamferedTriangle_PtsTest"

#procedure_name = "mcfCurve0"
#procedure_name = "mcfCurve1"

#directory = "../output/core_tests/"  # Adjust this path accordingly
directory = "../output"  # Adjust this path accordingly

fps = 10

png_time_steps = [] # a specified time step container for png export. If empty, a gif animation will be exported
svg_time_steps = [] # a specified time step container for svg export. If empty, all time steps will be exported

# a specified time step container for exporting a single png with opaque curve polylines (except the last one). If empty, no png will be exported
multi_png_time_steps = []
#multi_png_time_steps = [0, 10, 20, 30, 100, 200, 300, 400] 
#multi_png_time_steps = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100] 
#multi_png_time_steps = [0, 20, 40, 60, 80, 100] 
#multi_png_time_steps = [0, 10, 100] 
#multi_png_time_steps = [0, 4, 8, 12, 16, 20, 24, 26, 30] 

inner_curves_color_palette = [
    '#65107a', # dark purple
    #'#00ffff', # cyan
    #'#e3fcfc', # bright cyan
    '#1c7373', # dark teal
    #'#ff00ff', # magenta
    #'#32CD32', # lime green
]

inner_curves_line_styles = [
    (3, 1, 1, 1), # dash-dot
    (None, None),  # solid
    (None, None),  # solid
    (3, 1),       # dashed
    (1, 1),       # dotted
]

# =============================================================
#        Inscribed circle estimation
# -------------------------------------------------------------

inscribed_circle_parameters = {
    "armadillonewEvol_Pts2D3": {
        "center": (-3.0, 52.0),
        "radius": 20.0
    },
    "bunnynewEvol_Pts2D3": {
        "center": (-0.025, 0.08),
        "radius": 0.025
    },
    "maxPlancknewEvol_Pts2D3": {
        "center": (8.0, 85.0),
        "radius": 50.0
    },
    "nefertitinewEvol_Pts2D3": {
        #"center": (-75.0, -50.0),
        #"radius": 25.0
        "center": (-20.0, 90.0),
        "radius": 55.0
    }
}

inscribed_circle = None #inscribed_circle_parameters.get(procedure_name)

# =============================================================

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

background_img = None
legend_img = None

if os.path.exists(background_image_path) and os.path.exists(gdim2d_path):
    bbox_min, bbox_max, nx, ny, cell_size = read_gdim2d_file(gdim2d_path)
    background_img = imageio.imread(background_image_path)
    extent = [bbox_min[0], bbox_max[0], bbox_max[1], bbox_min[1]]

if os.path.exists(legend_image_path):
    legend_img = imageio.imread(legend_image_path)
    bbox_size = np.array(bbox_max) - np.array(bbox_min)
    clip_factor = 0.05
    legend_width = legend_img.shape[1]  # Get the image width in pixels
    clipped_width = int(legend_width * clip_factor)
    legend_img = legend_img[:, clipped_width:] # clip legend

    # Get the legend's aspect ratio (width/height)
    legend_height, legend_width = legend_img.shape[:2]
    legend_aspect_ratio = legend_width / legend_height

    # Define the height of the legend in terms of bbox size
    legend_height_in_bbox = 0.5 * bbox_size[1]
    legend_bbox_y_offset = (bbox_size[1] - legend_height_in_bbox) / 2.0
    legend_width_in_bbox = legend_height_in_bbox * legend_aspect_ratio  # Preserve aspect ratio

    extent_legend = [bbox_max[0],
                     bbox_max[0] + legend_width_in_bbox,
                     bbox_min[1] + legend_bbox_y_offset,
                     bbox_max[1] - legend_bbox_y_offset]


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
        # The files take form of "{procedure_name}_Inner{curve_id}_Evol_{time_id}.ply"
        # Extract the curve ID and time ID from the filename
        parts = f.split('_')

        # find part containing "Inner" to extract curve_id
        curve_id_part = [part for part in parts if "Inner" in part][0]
        curve_id = curve_id_part.replace('Inner', '')  # Extract curve ID (e.g., "0" from "Inner0")

        time_id = int(parts[-1].replace('Evol', '').replace('.ply', ''))  # Extract time ID
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
    # Check if outer or inner curves exist
    if outer_curves:
        # Use the first outer curve to set plot limits
        first_outer_curve = outer_curves[0]

        x_bounds = [np.min(first_outer_curve[:, 0]), np.max(first_outer_curve[:, 0])]
        y_bounds = [np.min(first_outer_curve[:, 1]), np.max(first_outer_curve[:, 1])]

        # expand the bounds by 10% to avoid clipping the curves
        y_range = y_bounds[1] - y_bounds[0]
        x_range = x_bounds[1] - x_bounds[0]
        y_bounds[0] -= 0.1 * y_range
        y_bounds[1] += 0.1 * y_range
        x_bounds[0] -= 0.1 * x_range
        x_bounds[1] += 0.1 * x_range

        ax_main.set_xlim(x_bounds)
        ax_main.set_ylim(y_bounds)
        
        ax_main.invert_yaxis()  # Reverse the direction of the y-axis
    elif inner_curves:
        # Use the first inner curve if outer curves are not available
        first_inner_curve_group = next(iter(inner_curves.values()))[0]
        ax_main.set_xlim(np.min(first_inner_curve_group[:, 0]), np.max(first_inner_curve_group[:, 0]))
        ax_main.set_ylim(np.min(first_inner_curve_group[:, 1]), np.max(first_inner_curve_group[:, 1]))
        ax_main.invert_yaxis()  # Reverse the direction of the y-axis
    else:
        raise ValueError("No curves (outer or inner) are available to set plot bounds.")

# Show legend
if legend_img is not None:
    ax_legend.imshow(legend_img, extent=extent_legend)

    # Load and overlay the symbol above the legend within the same axis
    legend_symbol = "./dGamma_200dpi.png"  # Path to your symbol file
    symbol_img = mpimg.imread(legend_symbol)
    
    # Get the symbol's aspect ratio (width/height)
    symbol_height, symbol_width = symbol_img.shape[:2]
    symbol_aspect_ratio = symbol_width / symbol_height

    # Define the height of the symbol in terms of bbox size
    symbol_height_in_bbox = 0.05 * bbox_size[1]  # % of the bbox height
    symbol_width_in_bbox = symbol_height_in_bbox * symbol_aspect_ratio  # Preserve aspect ratio
    extent_legend_size_x = extent_legend[1] - extent_legend[0]
    
    # Overlay the symbol in the upper portion of the same axis, preserving aspect ratio
    symbol_extent = [extent_legend[0] + 0.1 * extent_legend_size_x,  # Horizontally position the symbol
                     extent_legend[0] + 0.1 * extent_legend_size_x + symbol_width_in_bbox,  # Set the width based on aspect ratio
                     extent_legend[3] + 0.0 * symbol_height_in_bbox,  # Adjust vertical position (start)
                     extent_legend[3] + 1.0 * symbol_height_in_bbox  # Adjust vertical position (end)
                     ]
    
    # Display the symbol and legend with adjusted extent
    ax_legend.imshow(legend_img, extent=extent_legend, aspect='equal')
    ax_legend.imshow(symbol_img, extent=symbol_extent, aspect='equal')
    
    ax_legend.set_xlim(extent_legend[:2])
    ax_legend.set_ylim([extent_legend[2], extent_legend[3] + 0.1 * bbox_size[1]])  # Adjust ylim for the symbol placement

# =============================================

ax_main.set_aspect('equal')

# Initialize the lines for inner (magenta) and outer (black) curves
outer_line, = ax_main.plot([], [], 'k-', linewidth=2)  # Black outer curve
inner_line_handles = [ax_main.plot([], [], color='#65107a', linestyle='-', linewidth=2)[0] for _ in inner_curves]  # Different lines for each inner curve


# Initialize the line data
def init():
    outer_line.set_data([], [])
    for inner_lines in inner_line_handles:
        inner_lines.set_data([], [])
    return [outer_line] + inner_line_handles

# Update function for animation
def update(frame):
    # Update outer curve
    if outer_curves and outer_curves[frame].size > 0:
        x_outer, y_outer = outer_curves[frame].T
        outer_line.set_data(x_outer, y_outer)
    
    # Update inner curves
    for inner_lines, inner_curve_group in zip(inner_line_handles, inner_curves.values()):
        if inner_curve_group[frame].size > 0:
            x_inner, y_inner = inner_curve_group[frame].T
            inner_lines.set_data(x_inner, y_inner)

    if inscribed_circle:
        center = inscribed_circle["center"]
        radius = inscribed_circle["radius"]
        inscribed_circle_patch = patches.Circle(center, radius, edgecolor='purple', facecolor='none', linewidth=2)
        ax_main.add_patch(inscribed_circle_patch)
    
    # Update the title with the current time step using time_id
    ax_main.set_title(f"Time Step: {time_ids[frame]}", fontsize=14)

    return [outer_line] + inner_line_handles

# Update curves only for svg export (do not add them to ax_main)
def update_curves_only(frame):
    # Update outer curve
    if outer_curves and outer_curves[frame].size > 0:
        x_outer, y_outer = outer_curves[frame].T
        outer_line.set_data(x_outer, y_outer)
    
    # Update inner curves
    for inner_lines, inner_curve_group in zip(inner_line_handles, inner_curves.values()):
        if inner_curve_group[frame].size > 0:
            x_inner, y_inner = inner_curve_group[frame].T
            inner_lines.set_data(x_inner, y_inner)

    return [outer_line] + inner_line_handles

def plot_curves_with_increasing_opacity(frames):
    # Draw the base image and overlay the curves with opacity from 0 to 1 for the specified frames
    n_frames = len(frames)
    min_opacity = 0.2  # Minimum opacity for the curves

    outer_label = "F"  # Label for the outer curve
    inner_labels = [f"$G_{{{idx}}}$" for idx in range(1, len(inner_curves) + 1)]  # Labels for the inner curves (indexed from 1 in the article)

    for frame_idx, frame in enumerate(frames):

        alpha = min_opacity + (1.0 - min_opacity) * (frame_idx / (n_frames - 1))  # linear interpolation for opacity
        # alpha = min_opacity + (1.0 - min_opacity) * (1.0 - np.exp(-frame_idx / (n_frames - 1)))  # exponential interpolation for opacity

        # Update outer curve
        if outer_curves and outer_curves[frame].size > 0:
            x_outer, y_outer = outer_curves[frame].T
            outer_line.set_data(x_outer, y_outer)
            if frame_idx == n_frames - 1:  # Add label only for the last frame
                ax_main.plot(x_outer, y_outer, 'k-', linewidth=2, alpha=1, label=outer_label)[0]
            else:
                ax_main.plot(x_outer, y_outer, 'k-', linewidth=2, alpha=alpha)[0]

        # Update inner curves
        for idx, inner_curve_group in enumerate(inner_curves.values()):
            if inner_curve_group[frame].size > 0:
                x_inner, y_inner = inner_curve_group[frame].T
                inner_color = inner_curves_color_palette[idx]  # Use idx to access color palette
                inner_line_type = inner_curves_line_styles[idx]  # Use idx to access line style palette
                if frame_idx == n_frames - 1:  # Add labels only for the last frame
                    ax_main.plot(x_inner, y_inner, color=inner_color, linewidth=2, alpha=1, label=inner_labels[idx], dashes=inner_line_type)[0]
                else:
                    ax_main.plot(x_inner, y_inner, color=inner_color, linewidth=2, alpha=alpha, dashes=inner_line_type)[0]

    # Set title and add the legend on the last frame
    ax_main.set_title(f"{procedure_name}", fontsize=14)
    if frame_idx == n_frames - 1:
        ax_main.legend(loc='upper right')  # Position the legend in the upper right corner



# =============================================================================================
#                                  Visualization
# ---------------------------------------------------------------------------------------------

if png_time_steps or svg_time_steps or multi_png_time_steps:
    # Specialized image exports
    if png_time_steps:
        # Manually call the update function for each specified time step and export as PNG
        for frame_idx, time_step in enumerate(time_ids):
            if time_step in png_time_steps:
                update(frame_idx)  # Call update for this frame
                output_png_path = os.path.join(directory, f"{procedure_name}_Step{time_step}.png")
                print("Saving ", output_png_path)
                plt.savefig(output_png_path, dpi=300)  # Save the PNG for the current frame
    elif svg_time_steps:
        # Manually call the update function for each specified time step and export as SVG
        for frame_idx, time_step in enumerate(time_ids):
            if time_step in svg_time_steps:
                update_curves_only(frame_idx)  # Call update for this frame
                output_svg_path = os.path.join(directory, f"{procedure_name}_Step{time_step}.svg")
                print("Saving ", output_svg_path)
                
                # save outer_line and inner_line_handles to svg as paths in correct coordinates with flipped y-axis
                with open(output_svg_path, 'w') as f:
                    f.write(f'<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="100%" height="100%">\n')
                    f.write(f'<g transform="scale(1,-1) translate(0,-{extent[3]})">\n')
                    f.write(f'<path d="M')
                    for x, y in zip(outer_line.get_xdata(), outer_line.get_ydata()):
                        f.write(f'{x},{y} ')
                    f.write(f'Z" stroke="black" fill="none" stroke-width="2"/>\n')
                    for inner_lines in inner_line_handles:
                        f.write(f'<path d="M')
                        for x, y in zip(inner_lines.get_xdata(), inner_lines.get_ydata()):
                            f.write(f'{x},{y} ')
                        f.write(f'Z" stroke="#65107a" fill="none" stroke-width="2"/>\n')
                    f.write(f'</g>\n')
                    f.write(f'</svg>\n')
    elif multi_png_time_steps:
        # Draw the base image and overlay the curves with increasing opacity for the specified frames
        plot_curves_with_increasing_opacity(multi_png_time_steps)
        output_png_path = os.path.join(directory, f"{procedure_name}_CurvesOpacity.png")
        print("Saving ", output_png_path)
        plt.savefig(output_png_path, dpi=300)
else:
    # Create the animation
    # Determine the number of frames based on available data
    if outer_curves and len(outer_curves) > 0:
        n_frames = len(outer_curves)
    elif inner_curves and len(next(iter(inner_curves.values()))) > 0:
        # Check the first inner curve group for frames
        n_frames = len(next(iter(inner_curves.values())))
    else:
        raise ValueError("No curves (outer or inner) are available for animation.")
    ani = FuncAnimation(fig, update, frames=n_frames, init_func=init, blit=False, repeat=True)

    # Save the animation as a GIF
    output_gif_path = os.path.join(directory, f"{procedure_name}_animation.gif")
    ani.save(output_gif_path, writer='pillow', fps=fps)

# Show the animation
plt.show()
