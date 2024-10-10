import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import imageio
import json
import matplotlib.colors as mcolors
from matplotlib import cm
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
import numpy as np

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

# Function to parse the particle swarm log file
def parse_particle_swarm_log(json_file_path):
    with open(json_file_path, 'r') as f:
        data = json.load(f)
    
    # Convert point strings to tuples of coordinates
    for step, points in data.items():
        for key, value in points.items():
            x, y = map(int, value.strip("()").split(","))
            points[key] = (x, y)
    
    return data

# Helper function to convert grid indices to actual coordinates
def index_to_coordinates(coord, bbox_min, cell_size):
    return bbox_min[0] + coord[0] * cell_size, bbox_min[1] + coord[1] * cell_size

# Helper function to linearly interpolate between two colors
def interpolate_color(color1, color2, interp):
    c1_rgb = mcolors.hex2color(color1)
    c2_rgb = mcolors.hex2color(color2)
    interpolated_rgb = [(1 - interp) * c1 + interp * c2 for c1, c2 in zip(c1_rgb, c2_rgb)]
    return mcolors.to_hex(interpolated_rgb)

start_color = '#FF0000'  # Red
#start_color = '#000000'  # Black
#end_color = '#37a62b'    # Green
end_color = '#1cb0ff'    # Blue
#end_color = '#FF0000'    # Red

#edge_color = '#9e9e9e'  # Gray
edge_color = '#3b3b3b'  # Dark Gray

alpha_min, alpha_max = 1.0, 1.0  # Alpha scaling range

# Function to plot particle swarm paths with dynamic alpha and color scaling, and connect particles with line segments
def plot_particle_swarm(ax, swarm_data, bbox_min, cell_size):
    total_steps = len(swarm_data)
    previous_points = None  # To store the previous step's particle positions

    # First, draw the connecting lines behind the circles
    for step_index, (step, points) in enumerate(swarm_data.items()):
        # Calculate alpha for this step
        alpha = alpha_min + (alpha_max - alpha_min) * (step_index / (total_steps - 1))

        # Convert particle positions to actual coordinates and store for line drawing
        current_points = {}  # Store the current step's particle positions
        for key, (x, y) in points.items():
            coord_x, coord_y = index_to_coordinates((x, y), bbox_min, cell_size)
            current_points[key] = (coord_x, coord_y)

            # If we have previous points, connect them with a line
            if previous_points is not None and key in previous_points:
                prev_x, prev_y = previous_points[key]
                # Plot line connecting current and previous position of the same particle
                ax.plot([prev_x, coord_x], [prev_y, coord_y], color=edge_color, linewidth=1.5, alpha=alpha, zorder=1)

        # Update the previous_points with the current step's positions
        previous_points = current_points

    # Then, draw the circles on top of the lines
    for step_index, (step, points) in enumerate(swarm_data.items()):
        # Calculate alpha and color for this step
        alpha = alpha_min + (alpha_max - alpha_min) * (step_index / (total_steps - 1))
        color = interpolate_color(start_color, end_color, step_index / (total_steps - 1))

        # Convert particle positions to actual coordinates and draw particles
        for key, (x, y) in points.items():
            coord_x, coord_y = index_to_coordinates((x, y), bbox_min, cell_size)
            
            # Plot circle with white outline and interpolated color
            circle = patches.Circle(
                (coord_x, coord_y), 
                radius=0.3 * cell_size, 
                fill=True, 
                edgecolor=edge_color, 
                facecolor=color, 
                alpha=alpha,
                zorder=2  # zorder for circles
            )
            ax.add_patch(circle)

# Function to add a vertical color legend with discrete color values for each time step
def add_color_legend(fig, total_steps, start_color, end_color, alpha_min, alpha_max):
    # Create a list of RGBA colors that interpolate between start_color and end_color
    colors_with_alpha = [
        mcolors.to_rgba(interpolate_color(start_color, end_color, step / (total_steps - 1)), alpha_min + (alpha_max - alpha_min) * (step / (total_steps - 1)))
        for step in range(total_steps)
    ]
    
    # Create a discrete color map using the list of colors
    cmap = mcolors.ListedColormap(colors_with_alpha)
    
    # Normalization of time steps
    norm = Normalize(vmin=0, vmax=total_steps-1)

    # Create the color bar with discrete colors
    cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.7])  # Position for the vertical color bar
    colorbar = ColorbarBase(cbar_ax, cmap=cmap, norm=norm, orientation='vertical', boundaries=np.arange(-0.5, total_steps + 0.5, 1), spacing='uniform')

    # Add label for the color bar
    colorbar.set_label('Time Steps')

    # Set tick marks and labels to show each individual time step
    colorbar.set_ticks(np.arange(0, total_steps))
    colorbar.set_ticklabels([f'step {i}' for i in range(total_steps)])

# Function to draw white grid lines
def draw_grid_lines(ax, bbox_min, bbox_max, nx, ny, cell_size):
    # Draw vertical lines
    for i in range(nx):
        x = bbox_min[0] + i * cell_size
        ax.plot([x, x], [bbox_min[1], bbox_max[1]], color='white', linewidth=2, alpha=0.2)
    
    # Draw horizontal lines
    for i in range(ny):
        y = bbox_min[1] + i * cell_size
        ax.plot([bbox_min[0], bbox_max[0]], [y, y], color='white', linewidth=2, alpha=0.2)


# Main function to visualize the particle swarm paths
def visualize_particle_swarm(procedure_name, data_dir, grid_overlay=False):
    background_image_path = os.path.join(data_dir, f"{procedure_name}.png")
    gdim2d_path = os.path.join(data_dir, f"{procedure_name}.gdim2d")
    log_file_path = os.path.join(data_dir, f"{procedure_name}_ParticleSwarmLog.json")

    # Check if files exist
    if not (os.path.exists(background_image_path) and os.path.exists(gdim2d_path) and os.path.exists(log_file_path)):        
        raise FileNotFoundError("Required files are missing!")

    # Read gdim2d file
    bbox_min, bbox_max, nx, ny, cell_size = read_gdim2d_file(gdim2d_path)

    # Load background image
    background_img = imageio.imread(background_image_path)
    background_img = np.flipud(background_img)  # Flip the image data vertically

    # Adjust the extent using the real-world coordinates, scaled by pixel dimensions
    # The extent is calculated as bbox_min to (bbox_min + (number of cells * cell size))
    extent = [
        bbox_min[0],  # X-min
        bbox_min[0] + nx * cell_size,  # X-max
        bbox_min[1],  # Y-min
        bbox_min[1] + ny * cell_size  # Y-max
    ]

    # Parse particle swarm log file
    particle_swarm_data = parse_particle_swarm_log(log_file_path)

    # Plot the image and particle swarm paths
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(background_img, extent=extent, origin='upper')

    if grid_overlay:
        # Draw white grid lines based on grid dimensions
        draw_grid_lines(ax, bbox_min, bbox_max, nx, ny, cell_size)

    # Plot particle paths on top of the background image with alpha and color scaling
    plot_particle_swarm(ax, particle_swarm_data, bbox_min, cell_size)
    plt.xlabel("X")
    plt.ylabel("Y")

    # Add the color legend for time steps and alpha
    add_color_legend(fig, len(particle_swarm_data), start_color=start_color, end_color=end_color, alpha_min=alpha_min, alpha_max=alpha_max)

    #plt.title(f"Particle Swarm Paths: {procedure_name}")
    plt.show()

# Call the function with procedure name
visualize_particle_swarm("incompleteCircleDF", "../output", True)
