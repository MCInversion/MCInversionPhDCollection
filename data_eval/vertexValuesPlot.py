import json
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# Configuration
procedure_name = "equilibriumConcavePair16"
directory = "../output"  # Adjust this path accordingly
json_file = f"{directory}/{procedure_name}_log.json"

use_interactive_plot = False  # Toggle between slider and opacity visualization

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


# Opacity-based interpolation plot
def visualize_opacity_interpolation(min_opacity=0.2):
    selected_time_steps = [0, 2, 5, 20, 80, 200, 400, 500]  # Specific time steps to include
    n_steps = len(selected_time_steps)
    fig, ax = plt.subplots()
    ax.set_xlabel("Normalized Vertex Indices [0, 1]")
    ax.set_ylabel("Values")
    ax.set_title(f"{procedure_name} - Increasing Opacity")

    frame_idx = 0  # Initialize frame index
    for time_step in selected_time_steps:
        step_key = f"TimeStep_{time_step}"  # Adjusting to match the key format
        if step_key not in data:
            print(f"Warning: {step_key} not found in the data. Skipping.")
            continue

        step_data = data[step_key]
        alpha = min_opacity + (1.0 - min_opacity) * (frame_idx / (n_steps - 1)) if n_steps > 1 else 1.0

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
                        alpha=alpha,
                        label=f"{manifold_name} - {value_type}" if frame_idx == n_steps - 1 else None,
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
                        alpha=alpha,
                        label=f"{manifold_name} - {value_type}" if frame_idx == n_steps - 1 else None,
                    )
        
        frame_idx += 1  # Increment frame index

    # Add legend on the final frame
    if frame_idx == n_steps:
        ax.legend(loc="upper right")
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
