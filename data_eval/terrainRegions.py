import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

# Define boundaries for each example
boundaries = [
    # Example 1: Square with side 100 units
    [
        (0, 0),
        (100, 0),
        (100, 100),
        (0, 100)
    ],
    # Example 2: Triangle of roughly the same size as the 100-unit side square
    [
        (0, 0),
        (100, 0),
        (50, 86.6)  # Height of an equilateral triangle with side length 100 units
    ],
    # Example 3: More complicated concave polygon with edge length at least 100 units
    [
        (0, 0),
        (100, 0),
        (100, 50),
        (50, 50),
        (50, 100),
        (0, 100)
    ]
]

# Plot each boundary
fig, axs = plt.subplots(1, 3, figsize=(18, 6))

for i, boundary in enumerate(boundaries):
    polygon = Polygon(boundary, closed=True, edgecolor='black', facecolor='lightblue')
    axs[i].add_patch(polygon)
    axs[i].set_xlim(-10, 110)
    axs[i].set_ylim(-10, 110)
    axs[i].set_aspect('equal')
    axs[i].set_title(f'Example {i+1}')
    axs[i].grid(True)

plt.tight_layout()
plt.show()
