import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import imageio
import json
import random
import matplotlib.colors as mcolors
import numpy as np

# Function to convert hex color to RGB tuple
def hex_to_rgb(hex_color):
    return mcolors.hex2color(hex_color)

# Function to convert RGB tuple to hex color
def rgb_to_hex(rgb_color):
    return mcolors.to_hex(rgb_color)

# Function to linearly interpolate between two colors
def interpolate_color(base_color, end_color, interp):
    base_rgb = hex_to_rgb(base_color)
    end_rgb = hex_to_rgb(end_color)
    
    interpolated_rgb = [
        (1 - interp) * base_val + interp * end_val
        for base_val, end_val in zip(base_rgb, end_rgb)
    ]
    
    return rgb_to_hex(interpolated_rgb)

# Function to randomize color between base_color and end_color
def randomize_color_between(base_color, end_color):
    base_rgb = hex_to_rgb(base_color)
    end_rgb = hex_to_rgb(end_color)

    # Generate a random weight for interpolation
    random_weight = random.uniform(0, 1)

    # Interpolate between base_rgb and end_rgb
    randomized_rgb = [
        (1 - random_weight) * base_val + random_weight * end_val
        for base_val, end_val in zip(base_rgb, end_rgb)
    ]
    
    return rgb_to_hex(randomized_rgb)

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

# Recursive function to process each region in the JSON structure
def process_region(region, depth=0):
    processed_region = {
        "depth": depth,
        "bounds": (region['minX'], region['maxX'], region['minY'], region['maxY']),
        "is_maximum": region.get("Maximum found!", False),  # Handle maximum flag
        "children": []
    }
    
    # Process each child region recursively
    for child in region.get('children', []):
        processed_region['children'].append(process_region(child, depth + 1))
    
    return processed_region

# Function to parse the JSON log file
def parse_log_file(json_file_path):
    with open(json_file_path, 'r') as f:
        data = json.load(f)  # Load the JSON data (ensure the content is JSON formatted correctly)
    
    # Ensure data is a list or dictionary (depending on your actual JSON structure)
    if isinstance(data, dict):
        return [process_region(data)]
    elif isinstance(data, list):
        return [process_region(region) for region in data]
    else:
        raise ValueError("Unexpected JSON format.")


# Helper function to convert grid indices to coordinate bounds (with flipped Y-axis)
def index_to_coord(min_idx, max_idx, bbox_min, bbox_max, cell_size, is_y=False):
    # Calculate the minimum and maximum coordinates
    min_coord = bbox_min + (min_idx - 1) * cell_size #+ cell_size * (1 if is_y else -1)
    max_coord = bbox_min + (max_idx - 1) * cell_size #+ cell_size * (2 if is_y else -1)
    # # Flip the Y-axis if needed
    # if is_y:
    #     max_coord = bbox_max - min_idx * cell_size
    #     min_coord = bbox_max - max_idx * cell_size

    return min_coord, max_coord

# Function to calculate the maximum depth of regions
def calculate_max_depth(regions):
    if not regions:
        return 0
    return max(region["depth"] + calculate_max_depth(region["children"]) for region in regions)

# Helper function to flatten the hierarchical regions into a list with depths
def flatten_regions(regions):
    flat_regions = []

    def recurse(region_list, current_depth):
        for region in region_list:
            flat_regions.append((current_depth, region))  # Add region with its depth
            if region["children"]:
                recurse(region["children"], current_depth + 1)

    recurse(regions, 0)  # Start with depth 0
    return flat_regions

# Plot interpolation settings
alpha_base = 0.05
alpha_max = 0.8
width_base = 1
width_max = 1

base_color = '#255c27'
end_color = '#511d6b'

# Function to plot regions from lowest to highest depth, with optional min_depth
def plot_hierarchical_regions(ax, regions, bbox_min, bbox_max, cell_size, max_depth, min_depth=0, plot_maxima=False):
    print("min_depth: ", min_depth, ", max_depth: ", max_depth)
    use_single_depth = max_depth == min_depth

    # Flatten the regions into a list with their depths
    flat_regions = flatten_regions(regions)
    
    # Sort regions by depth (from lowest to highest)
    flat_regions.sort(key=lambda x: x[0])

    # Plot each region in the order of increasing depth
    for depth, region in flat_regions:
        if depth < min_depth:
            continue  # Skip regions shallower than the specified min_depth

        if use_single_depth and depth != min_depth:
            continue  # Skip regions that are not exactly at the specified min_depth

        if region["bounds"]:
            minX, maxX, minY, maxY = region["bounds"]
            min_coord_x, max_coord_x = index_to_coord(minX, maxX, bbox_min[0], bbox_max[0], cell_size, False)
            min_coord_y, max_coord_y = index_to_coord(minY, maxY, bbox_min[1], bbox_max[1], cell_size, True)

            interp = 1 if max_depth == min_depth else 1 - (depth - min_depth) / (max_depth - min_depth)

            alpha = ((1 - interp) * alpha_base + interp * alpha_max) * 1
            width = 2
            
            # Interpolate color between base_color and end_color
            #color = interpolate_color(base_color, end_color, interp)
            color = randomize_color_between(base_color, end_color)
            
            # Use red color for regions marked as maximum, black otherwise
            color = 'red' if plot_maxima and region["is_maximum"] else color
            alpha = 1 if plot_maxima and region["is_maximum"] else alpha

            rect = patches.Rectangle((min_coord_x, min_coord_y), max_coord_x - min_coord_x,
                                     max_coord_y - min_coord_y, linewidth=width, edgecolor=color, 
                                     facecolor='none', alpha=alpha)
            ax.add_patch(rect)

# Function to plot maximal regions from min_depth to max_depth
def plot_maximal_regions(ax, regions, bbox_min, bbox_max, cell_size, max_depth, min_depth=0):
    print("min_depth: ", min_depth, ", max_depth: ", max_depth)
    use_single_depth = max_depth == min_depth

    # Flatten the regions into a list with their depths
    flat_regions = flatten_regions(regions)
    
    # Sort regions by depth (from lowest to highest)
    flat_regions.sort(key=lambda x: x[0])

    # Plot each region in the order of increasing depth
    for depth, region in flat_regions:
        if depth < min_depth:
            continue  # Skip regions shallower than the specified min_depth

        if use_single_depth and depth != min_depth:
            continue  # Skip regions that are not exactly at the specified min_depth

        if region["bounds"] and region["is_maximum"]:
            minX, maxX, minY, maxY = region["bounds"]
            min_coord_x, max_coord_x = index_to_coord(minX, maxX, bbox_min[0], bbox_max[0], cell_size, False)
            min_coord_y, max_coord_y = index_to_coord(minY, maxY, bbox_min[1], bbox_max[1], cell_size, True)
            #center_x, center_y = 0.5 * (min_coord_x + max_coord_x), 0.5 * (min_coord_y + max_coord_y)
            rect = patches.Rectangle((min_coord_x, min_coord_y), max_coord_x - min_coord_x,
                                     max_coord_y - min_coord_y, linewidth=2, edgecolor='red', 
                                     facecolor='none', alpha=1)
            ax.add_patch(rect)
            #center = patches.Circle(xy=(center_x, center_y), radius=0.25*cell_size, facecolor='none', edgecolor='red')
            #ax.add_patch(center)


# Function to print the depths of all regions
def print_region_depths(regions):
    # Helper function to recursively print the depth of regions
    def recurse(region_list, current_depth):
        for region in region_list:
            print(f"Depth: {current_depth}, Bounds: {region.get('bounds', None)}")
            if region["children"]:
                recurse(region["children"], current_depth + 1)

    recurse(regions, 0)  # Start with depth 0

grid_color = '#ffffff'  # White
#grid_color = '#3b3b3b'  # Dark Gray

# Function to draw white grid lines
def draw_grid_lines(ax, bbox_min, bbox_max, nx, ny, cell_size):
    # Draw vertical lines
    for i in range(nx):
        x = bbox_min[0] + i * cell_size
        ax.plot([x, x], [bbox_min[1], bbox_max[1]], color=grid_color, linewidth=2, alpha=0.2)
    
    # Draw horizontal lines
    for i in range(ny):
        y = bbox_min[1] + i * cell_size
        ax.plot([bbox_min[0], bbox_max[0]], [y, y], color=grid_color, linewidth=2, alpha=0.2)

# Main function to visualize the progress
def visualize_inscribed_circle_calculation(procedure_name, data_dir, target_depth=None, grid_overlay=False):
    background_image_path = os.path.join(data_dir, f"{procedure_name}.png")
    gdim2d_path = os.path.join(data_dir, f"{procedure_name}.gdim2d")
    log_file_path = os.path.join(data_dir, f"{procedure_name}_HierarchicalLog.json")

    # Check if files exist
    if not (os.path.exists(background_image_path) and os.path.exists(gdim2d_path) and os.path.exists(log_file_path)):        
        print("background_image_path: ", background_image_path)
        print("gdim2d_path: ", gdim2d_path)
        print("log_file_path: ", log_file_path)
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

    # Parse hierarchical log file
    hierarchical_regions = parse_log_file(log_file_path)

    # Plot the image and hierarchical regions
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(background_img, extent=extent, origin='upper')

    if grid_overlay:
        # Draw white grid lines based on grid dimensions
        draw_grid_lines(ax, bbox_min, bbox_max, nx, ny, cell_size)

    # Calculate the maximum depth of the hierarchy
    max_depth = calculate_max_depth(hierarchical_regions)
    min_depth = 0
    if target_depth is not None:
        min_depth = target_depth
        max_depth = target_depth

    print_region_depths(hierarchical_regions)

    # Plot using the max depth
    #plot_hierarchical_regions(ax, hierarchical_regions, bbox_min, bbox_max, cell_size, max_depth=max_depth, min_depth=min_depth)
    plot_maximal_regions(ax, hierarchical_regions, bbox_min, bbox_max, cell_size, max_depth=max_depth, min_depth=min_depth)

    plt.title(f"Inscribed Circle Calculation Progress: {procedure_name}")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.show()

# Call the function with procedure name
visualize_inscribed_circle_calculation("incompleteCircleDF", "../output", target_depth=4, grid_overlay=True)
