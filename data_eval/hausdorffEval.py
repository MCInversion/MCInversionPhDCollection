import matplotlib.pyplot as plt

# Path to the file containing Hausdorff distance measurements
file_path = 'HausdorffDistEvolution_Results.txt'

# Define the scaling factors for the specific procedures
scaling_factors = {
    "bunnyRepulsive": 1/10,  # Scaling factor for "Repulsive" to "Base"
    "bunnyDanielLSW": 1/1000,  # Scaling factor for "Daniel" to "Base"
    "bunnyLSWObstacle": 1,  # Base
    "bunnyLSWFullWrap": 1  # Base
}

# Initialize data structures for storing parsed data
procedure_data = {}

# Read and parse the file
with open(file_path, 'r') as file:
    lines = file.readlines()
    current_procedure = ''
    for line in lines:
        line = line.strip()
        if line.startswith('Procedure:'):
            current_procedure = line.split('Procedure:')[1].split('Hausdorff')[0].strip()
            procedure_data[current_procedure] = []
        elif line.startswith('step'):
            _, distance = line.split(': ')
            # Scale the distance if the current procedure has a defined scaling factor
            scaled_distance = float(distance) * scaling_factors.get(current_procedure, 1)
            procedure_data[current_procedure].append(scaled_distance)

# Plotting the data
plt.figure(figsize=(8, 4))
for procedure, distances in procedure_data.items():
    # Plot line with line width of 2
    line, = plt.plot(distances, label=procedure, linewidth=2)
    # Mark start and end points with a circle of the same color as the line
    plt.scatter([0], [distances[0]], color=line.get_color(), s=40)  # Start point
    plt.scatter([len(distances)-1], [distances[-1]], color=line.get_color(), s=40)  # End point

plt.title('Evolution of Mesh to Point Cloud Hausdorff Distance')
plt.xlabel('Time Step')
plt.ylabel('Hausdorff Distance')
plt.legend()
plt.grid(True)
plt.show()
