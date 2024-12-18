import json
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# Configuration
procedure_name = "equilibriumConcavePair19"
#procedure_name = "equilibriumPair0"
#procedure_name = "equilibriumPair1"
#procedure_name = "equilibriumPair2"
#procedure_name = "equilibriumPair3"

directory = "../output"  # Adjust this path accordingly
json_file = f"{directory}/{procedure_name}_log.json"

use_interactive_plot = False  # Toggle between slider and opacity visualization
force_vertex_indices_on_x_axis = False  # Use vertex indices instead of arc lengths on x-axis

# Line style configuration
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
chosen_lwd = 2  # Line width for all curves

# Load JSON data
with open(json_file, "r") as f:
    data = json.load(f)

# Extract all time steps
time_steps = sorted(data.keys(), key=lambda x: int(x.split("_")[-1]))

# Normalize vertex indices
def normalize_vertex_indices(vertex_indices):
    min_idx = min(vertex_indices)
    max_idx = max(vertex_indices)
    return [(idx - min_idx) / (max_idx - min_idx) for idx in vertex_indices]

# Helper function to extract manifold index
def extract_manifold_idx(manifold_name):
    match = re.search(r'\d+', manifold_name)  # Extract first numeric group
    return int(match.group()) if match else 0

# Helper function to normalize values to [0, 1]
def normalize_values(values):
    min_val = min(values)
    max_val = max(values)
    return [(v - min_val) / (max_val - min_val) for v in values]

# Helper function to get vertex ordering based on a given array
def get_vertex_ordering(input_list, key):
    return [i for i, _ in sorted(enumerate(input_list), key=lambda x: x[1][key])]

# Helper function to check uniqueness of values with a tolerance
def check_uniqueness_with_tolerance(values, tolerance=1e-6):
    sorted_values = sorted(values)  # Sort the values to check neighbors efficiently

    for i in range(1, len(sorted_values)):
        if abs(sorted_values[i] - sorted_values[i - 1]) < tolerance:
            print(f"Duplicate values found: {sorted_values[i - 1]} and {sorted_values[i]} within tolerance {tolerance}")
            return False

    return True

# Slider-based interactive plot with fixed plot range
def interactive_plot():
    n_steps = len(time_steps)
    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.25)

    # Determine global min/max values across all time steps
    global_min, global_max = float('inf'), float('-inf')
    for step_key in time_steps:
        step_data = data[step_key]
        for manifold_name, values in step_data.items():
            for value_type, value_list in values.items():
                values = [v["value"] for v in value_list]
                global_min = min(global_min, min(values))
                global_max = max(global_max, max(values))

    # Plotting function for slider
    def plot_values(time_step_id):
        ax.clear()
        step_key = time_steps[time_step_id]
        step_data = data[step_key]

        for manifold_name, values in step_data.items():
            manifold_idx = extract_manifold_idx(manifold_name)
            for value_type, value_list in values.items():
                vertex_indices = [v["vertexIndex"] for v in value_list]
                values = [v["value"] for v in value_list]
                normalized_indices = normalize_vertex_indices(vertex_indices)

                if manifold_idx == 0:
                    # Outer curve (black line)
                    ax.plot(
                        normalized_indices,
                        values,
                        color='black',
                        linewidth=chosen_lwd,
                        label=f"{manifold_name} - {value_type}",
                    )
                else:
                    # Inner curves with styles and colors
                    curve_idx = (manifold_idx - 1) % len(inner_curves_color_palette)
                    ax.plot(
                        normalized_indices,
                        values,
                        color=inner_curves_color_palette[curve_idx],
                        linewidth=chosen_lwd,
                        dashes=inner_curves_line_styles[curve_idx],
                        label=f"{manifold_name} - {value_type}",
                    )

        ax.set_xlabel("Normalized Vertex Indices [0, 1]")
        ax.set_ylabel("Values")
        ax.set_ylim(global_min, global_max)  # Fixed plot range
        ax.legend()
        ax.set_title(f"Time Step: {time_step_id}")
        plt.draw()

    # Initial plot
    plot_values(0)

    # Slider setup
    ax_slider = plt.axes([0.2, 0.1, 0.65, 0.03])
    slider = Slider(ax_slider, "Time Step", 0, n_steps - 1, valinit=0, valstep=1)

    def update(val):
        plot_values(int(slider.val))

    slider.on_changed(update)
    plt.show()


def get_sorted_values(values, manifold_name):
    # Extract "arcLength" and compute the vertex order
    arc_lengths = values.get("arcLength", [])

    # Check if values contains other values besides "arcLength"                    
    next_not_arc_length = next((k for k in values.keys() if k != "arcLength"), None)
    if next_not_arc_length is None:
        print(f"ERROR: No other values found in {manifold_name}. Something went wrong!")
        return(None) 

    if arc_lengths and not force_vertex_indices_on_x_axis:
        # Working with arc lengths
        vertex_order = get_vertex_ordering(arc_lengths, "value")
        sorted_arc_lengths = [arc_lengths[i]["value"] for i in vertex_order]
        sorted_values = []
        # We need to extend the value lists by two values to extrapolate the first edge
        max_vertex_arc_length = sorted_arc_lengths[-1]
        sub_max_vertex_arc_length = sorted_arc_lengths[-2]
        extrapolated_total_arc_length = 2 * max_vertex_arc_length - sub_max_vertex_arc_length
        extrapolated_difference = extrapolated_total_arc_length - max_vertex_arc_length
        ref_edge_length = extrapolated_difference + sorted_arc_lengths[0]

        sorted_arc_lengths.insert(0, 0.0)
        sorted_arc_lengths.append(extrapolated_total_arc_length)

        # extended_vertex_order = vertex_order.copy()
        # extended_vertex_order.insert(0, -1)
        # extended_vertex_order.append(len(arc_lengths))

        arc_length_param = extrapolated_difference / ref_edge_length

        # Create a new sorted value list with respective value types (without arc lengths)
        for value_type, value_list in values.items():
            if value_type == "arcLength":
                continue
                        
            sorted_value_list = [value_list[i]["value"] for i in vertex_order]

            # what value should be at the arc length 0.0?
            prev_value = sorted_value_list[-1]
            next_value = sorted_value_list[0]
            last_edge_interpolated_value = prev_value + arc_length_param * (next_value - prev_value)
            sorted_value_list.insert(0, last_edge_interpolated_value)
            sorted_value_list.append(last_edge_interpolated_value)

            value_type_sorted_values_pair = (value_type, sorted_value_list)
            sorted_values.append(value_type_sorted_values_pair)

        normalized_x = normalize_values(sorted_arc_lengths)
        return normalized_x, sorted_values

    # Fallback to vertex indices stored in every value list
    vertex_order = get_vertex_ordering(next_not_arc_length, "vertexIndex")

    # Create a new sorted value list with respective value (without arc lengths)
    sorted_values = []
    for value_type, value_list in values.items():
        if value_type == "arcLength":
            continue

        sorted_value_list = [value_list[i]["value"] for i in vertex_order]
        value_type_sorted_values_pair = (value_type, sorted_value_list)
        sorted_values.append(value_type_sorted_values_pair)

    normalized_x = normalize_vertex_indices(range(len(vertex_order)))
    return normalized_x, sorted_values


# Opacity-based interpolation plot
def visualize_opacity_interpolation(min_opacity=0.2):
    selected_time_steps = [1, 2, 8, 20, 50, 100, 250, 500, 1000, 1500, 2000, 2500]  # Selected time steps (1370)
    n_steps = len(selected_time_steps)
    fig, ax = plt.subplots()
    ax.set_ylabel("Values")
    ax.set_title(f"{procedure_name} - Increasing Opacity")

    frame_idx = 0  # Initialize frame index
    for time_step in selected_time_steps:
        step_key = f"TimeStep_{time_step}"
        if step_key not in data:
            print(f"ERROR: {step_key} not found in the data. Something went wrong!")
            break

        step_data = data[step_key]
        alpha = min_opacity + (1.0 - min_opacity) * (frame_idx / (n_steps - 1)) if n_steps > 1 else 1.0

        for manifold_name, values in step_data.items():
            manifold_idx = extract_manifold_idx(manifold_name)

            normalized_x, sorted_value_list = get_sorted_values(values, manifold_name)

            # Plot all other value types using the determined vertex order
            for value_type, y_values in sorted_value_list:

                if value_type == "arcLength":
                    print(f"ERROR: {value_type} found in {manifold_name}. Something went wrong!")
                    break

                if manifold_idx == 0:
                    # Outer curve (black line)
                    ax.plot(
                        normalized_x,
                        y_values,
                        color='black',
                        linewidth=chosen_lwd,
                        alpha=alpha,
                        label=f"{manifold_name} - {value_type}" if frame_idx == n_steps - 1 else None,
                    )
                else:
                    # Inner curves with styles and colors
                    curve_idx = (manifold_idx - 1) % len(inner_curves_color_palette)
                    ax.plot(
                        normalized_x,
                        y_values,
                        color=inner_curves_color_palette[curve_idx],
                        linewidth=chosen_lwd,
                        dashes=inner_curves_line_styles[curve_idx],
                        alpha=alpha,
                        label=f"{manifold_name} - {value_type}" if frame_idx == n_steps - 1 else None,
                    )

        frame_idx += 1  # Increment frame index

    # Add legend on the final frame
    if frame_idx == n_steps:
        ax.legend(loc="upper right")

    # Save and show the plot
    output_png_path = f"{directory}/{procedure_name}_ValueListsOpacity.png"
    print("Saving ", output_png_path)
    plt.savefig(output_png_path, dpi=300)
    plt.show()

# Main execution
if __name__ == "__main__":
    if use_interactive_plot:
        interactive_plot()
    else:
        visualize_opacity_interpolation(min_opacity=0.2)
