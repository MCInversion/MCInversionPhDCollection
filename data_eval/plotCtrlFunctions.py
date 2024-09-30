import matplotlib.pyplot as plt

# Data for outer and inner epsilon and eta
OUTER_DEBUG_EPSILON = [
	0.00152452,
	0.00124351,
	0.00122836,
	0.00122923,
	0.00124828,
	0.0011335,
	0.00106063,
	0.00092239,
	0.000814231,
	0.000768413,
	0.00095462,
	0.000617901,
	0.000558107,
	0.000565831,
	0.000616678,
	0.000909009,
	0.000760327,
	0.000811513,
	0.000917136,
	0.00105352
]
OUTER_DEBUG_ETA = [
	-0.0416522,
	-0.0362566,
	-0.0357816,
	-0.0357535,
	-0.0365151,
	-0.0405557,
	-0.0340123,
	-0.0315753,
	-0.0295651,
	-0.0286857,
	-0.0329489,
	-0.0261125,
	-0.0241112,
	-0.024615,
	-0.0259752,
	-0.0321462,
	-0.0285401,
	-0.0295423,
	-0.0315112,
	-0.0337706
]
INNER_DEBUG_EPSILON = [
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0
]
INNER_DEBUG_ETA = [
	-0.0256608,
	-0.0226758,
	-0.0227508,
	-0.0220763,
	-0.0208097,
	-0.0221332,
	-0.0227605,
	-0.0255329,
	-0.0254327,
	-0.0297939,
	-0.0259161,
	-0.0270085,
	-0.0269555,
	-0.0291751,
	-0.0305258,
	-0.0306042,
	-0.0349974,
	-0.0342155,
	-0.0333877,
	-0.0332035,
	-0.0322579,
	-0.0331493,
	-0.0333697,
	-0.034414,
	-0.0364545,
	-0.0308066,
	-0.0305439,
	-0.0291661,
	-0.0269307,
	-0.0269896,
	-0.0259857,
	-0.0260134
]

# Create a figure with two subplots for epsilon and eta
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

# Plot epsilon values for outer and inner manifolds
ax1.plot(OUTER_DEBUG_EPSILON, label='Outer Epsilon', color='#45caff', marker='o')
ax1.plot(INNER_DEBUG_EPSILON, label='Inner Epsilon', color='blue', marker='o')
ax1.set_title('Epsilon Control Weight')
ax1.set_xlabel('Vertex Index')
ax1.set_ylabel('Epsilon')
# Set y-limits dynamically based on the data
min_y1 = min(min(OUTER_DEBUG_EPSILON), min(INNER_DEBUG_EPSILON))
max_y1 = max(max(OUTER_DEBUG_EPSILON), max(INNER_DEBUG_EPSILON))
ax1.set_ylim(min_y1, max_y1)
ax1.grid(True)
ax1.legend()

# Get extremal indices for epsilon
outer_min_idx1 = OUTER_DEBUG_EPSILON.index(min(OUTER_DEBUG_EPSILON))
outer_max_idx1 = OUTER_DEBUG_EPSILON.index(max(OUTER_DEBUG_EPSILON))
inner_min_idx1 = INNER_DEBUG_EPSILON.index(min(INNER_DEBUG_EPSILON))
inner_max_idx1 = INNER_DEBUG_EPSILON.index(max(INNER_DEBUG_EPSILON))

# Draw dashed vertical lines at extremal points for epsilon
ax1.axvline(x=outer_min_idx1, ymin=0, ymax=(min(OUTER_DEBUG_EPSILON)-min_y1)/(max_y1-min_y1), color='#45caff', linestyle='--')
ax1.axvline(x=outer_max_idx1, ymin=0, ymax=(max(OUTER_DEBUG_EPSILON)-min_y1)/(max_y1-min_y1), color='#45caff', linestyle='--')
ax1.axvline(x=inner_min_idx1, ymin=0, ymax=(min(INNER_DEBUG_EPSILON)-min_y1)/(max_y1-min_y1), color='blue', linestyle='--')
ax1.axvline(x=inner_max_idx1, ymin=0, ymax=(max(INNER_DEBUG_EPSILON)-min_y1)/(max_y1-min_y1), color='blue', linestyle='--')

# Plot eta values for outer and inner manifolds
ax2.plot(OUTER_DEBUG_ETA, label='Outer Eta', color='#45caff', marker='o')
ax2.plot(INNER_DEBUG_ETA, label='Inner Eta', color='blue', marker='o')
ax2.set_title('Eta Control Weight')
ax2.set_xlabel('Vertex Index')
ax2.set_ylabel('Eta')
# Set y-limits dynamically based on the data
min_y2 = min(min(OUTER_DEBUG_ETA), min(INNER_DEBUG_ETA))
max_y2 = max(max(OUTER_DEBUG_ETA), max(INNER_DEBUG_ETA))
ax2.set_ylim(min_y2, max_y2)
ax2.grid(True)
ax2.legend()

# Get extremal indices for eta
outer_min_idx2 = OUTER_DEBUG_ETA.index(min(OUTER_DEBUG_ETA))
outer_max_idx2 = OUTER_DEBUG_ETA.index(max(OUTER_DEBUG_ETA))
inner_min_idx2 = INNER_DEBUG_ETA.index(min(INNER_DEBUG_ETA))
inner_max_idx2 = INNER_DEBUG_ETA.index(max(INNER_DEBUG_ETA))

# Draw dashed vertical lines at extremal points for eta
ax2.axvline(x=outer_min_idx2, ymin=0, ymax=(min(OUTER_DEBUG_ETA)-min_y2)/(max_y2-min_y2), color='#45caff', linestyle='--')
ax2.axvline(x=outer_max_idx2, ymin=0, ymax=(max(OUTER_DEBUG_ETA)-min_y2)/(max_y2-min_y2), color='#45caff', linestyle='--')
ax2.axvline(x=inner_min_idx2, ymin=0, ymax=(min(INNER_DEBUG_ETA)-min_y2)/(max_y2-min_y2), color='blue', linestyle='--')
ax2.axvline(x=inner_max_idx2, ymin=0, ymax=(max(INNER_DEBUG_ETA)-min_y2)/(max_y2-min_y2), color='blue', linestyle='--')

# Adjust layout to prevent overlap
plt.tight_layout()

# Display the plot
plt.show()

