import matplotlib.pyplot as plt

# Data for outer and inner epsilon and eta
OUTER_DEBUG_EPSILON = [
	0.00117694,
	0.000962369,
	0.000946529,
	0.000925103,
	0.000912722,
	0.000852732,
	0.000814528,
	0.000713163,
	0.000631303,
	0.000600952,
	0.000600952,
	0.000477927,
	0.000436568,
	0.000430658,
	0.000430658,
	0.000430658,
	0.000430658,
	0.000430658,
	0.000430658,
	0.000430658
]
OUTER_DEBUG_ETA = [
	-0.0372932,
	-0.0323466,
	-0.0311269,
	-0.0309917,
	-0.0316242,
	-0.0352184,
	-0.0299727,
	-0.0275399,
	-0.0258678,
	-0.0257096,
	-0.0313531,
	-0.023183,
	-0.0212318,
	-0.0212569,
	-0.0259466,
	-0.0288871,
	-0.029104,
	-0.0265734,
	-0.0215035,
	-0.0132155
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
	0.383417,
	0.583591,
	0.660167,
	0.6526,
	0.709565,
	0.679528,
	0.707673,
	0.646911,
	0.462061,
	0.214128,
	-0.0700699,
	-0.359203,
	-0.633229,
	-0.856458,
	-1.02494,
	-1.13646,
	-1.18291,
	-1.17048,
	-1.12068,
	-1.07521,
	-1.0613,
	-1.07255,
	-1.1154,
	-1.16649,
	-1.18283,
	-1.13936,
	-1.03156,
	-0.864549,
	-0.649714,
	-0.402953,
	-0.142664,
	0.110982
]


# Create a figure with two subplots for epsilon and eta
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

# Plot epsilon values for outer and inner manifolds
ax1.plot(OUTER_DEBUG_EPSILON, label='Outer Epsilon', color='blue', marker='o')
ax1.plot(INNER_DEBUG_EPSILON, label='Inner Epsilon', color='green', marker='x')
ax1.set_title('Epsilon Control Weight')
ax1.set_xlabel('Vertex Index')
ax1.set_ylabel('Epsilon')
# Set y-limits dynamically based on the data
ax1.set_ylim(min(min(OUTER_DEBUG_EPSILON), min(INNER_DEBUG_EPSILON)), max(max(OUTER_DEBUG_EPSILON), max(INNER_DEBUG_EPSILON)))
ax1.grid(True)
ax1.legend()

# Plot eta values for outer and inner manifolds
ax2.plot(OUTER_DEBUG_ETA, label='Outer Eta', color='blue', marker='o')
ax2.plot(INNER_DEBUG_ETA, label='Inner Eta', color='green', marker='x')
ax2.set_title('Eta Control Weight')
ax2.set_xlabel('Vertex Index')
ax2.set_ylabel('Eta')
# Set y-limits dynamically based on the data
ax2.set_ylim(min(min(OUTER_DEBUG_ETA), min(INNER_DEBUG_ETA)), max(max(OUTER_DEBUG_ETA), max(INNER_DEBUG_ETA)))
ax2.grid(True)
ax2.legend()

# Adjust layout to prevent overlap
plt.tight_layout()

# Display the plot
plt.show()
