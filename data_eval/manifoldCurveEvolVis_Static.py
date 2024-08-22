import os
import matplotlib.pyplot as plt
import numpy as np
import imageio.v2 as imageio  # Use imageio.v2 to avoid the deprecation warning
from collections import defaultdict

# Define the procedure name and the path to your directory containing the files
procedure_name = "ShrinkWrappingACirclePointCloud_NoInnerCurveNoRemeshing"
directory = "../output/core_tests/"  # Adjust this path accordingly

# Load the background image (distance field PNG) if it exists
background_image_path = os.path.join(directory, f"{procedure_name}TargetDF.png")
background_img = None
if os.path.exists(background_image_path):
    background_img = imageio.imread(background_image_path)
else:
    print(f"Background image {background_image_path} not found. Continuing without background.")

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
inner_ply_files = defaultdict(list)

for f in sorted(os.listdir(directory)):
    if f.startswith(f"{procedure_name}_Inner_") and f.endswith(".ply"):
        inner_id = f.split('_')[0]
        inner_ply_files[inner_id].append(os.path.join(directory, f))

outer_ply_files = sorted([
    os.path.join(directory, f) for f in os.listdir(directory)
    if f.startswith(f"{procedure_name}_Outer_Evol") and f.endswith(".ply")
])

# Read all vertex positions from the PLY files for inner and outer curves
inner_curves = {inner_id: [read_ply(ply_file) for ply_file in files]
                for inner_id, files in inner_ply_files.items()}
outer_curves = [read_ply(ply_file) for ply_file in outer_ply_files]

# Set up the plot
fig, ax = plt.subplots()

# Display the background image if it was successfully loaded
if background_img is not None:
    ax.imshow(background_img, extent=[-1.5, 1.5, -1.5, 1.5])

# Set up the plot limits and labels
ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1.5, 1.5)
ax.set_aspect('equal')

# Define transparency levels, from 20% to 100%
alphas = np.linspace(0.2, 1.0, len(outer_curves))

# Plot the outer curves with increasing transparency
for i, outer_curve in enumerate(outer_curves):
    if outer_curve.size > 0:
        x_outer, y_outer = outer_curve.T
        ax.plot(x_outer, y_outer, 'k-', linewidth=1, alpha=alphas[i])

# Plot the inner curves with increasing transparency
for inner_id, inner_curve_group in inner_curves.items():
    for i, inner_curve in enumerate(inner_curve_group):
        if inner_curve.size > 0:
            x_inner, y_inner = inner_curve.T
            ax.plot(x_inner, y_inner, 'b-', linewidth=1, alpha=alphas[i], label=f'{inner_id} Step {i}')


# Save the final image
output_image_path = os.path.join(directory, f"{procedure_name}_overlay.png")
plt.savefig(output_image_path, dpi=300)

# Show the final plot
plt.show()
