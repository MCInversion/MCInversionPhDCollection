import numpy as np

# def read_ply_vertex(file_path):
#     vertices = []
#     try:
#         with open(file_path, 'r') as f:
#             lines = f.readlines()

#         num_vertices = 0

#         # Parse header to get number of vertices
#         for i, line in enumerate(lines):
#             if line.startswith('element vertex'):
#                 num_vertices = int(line.split()[2])
#             elif line.startswith('end_header'):
#                 header_end = i
#                 break

#         # Read vertices
#         for line in lines[header_end + 1:header_end + 1 + num_vertices]:
#             coords = list(map(float, line.strip().split()))
#             vertices.append((coords[0], coords[1]))

#         # Ensure the polyline is closed by adding the starting point at the end
#         if vertices[-1] != vertices[0]:
#             vertices.append(vertices[0])

#     except Exception as e:
#         print(f"Error reading {file_path}: {e}")
#     return np.array(vertices)

def read_ply_edge(file_path):
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


import matplotlib.pyplot as plt

# Set the file path
file_path = "../output/bunnynewEvol_Pts2D3_Outer_Evol_80.ply"

# Read the polyline using the read_ply function
#polyline1 = read_ply_vertex(file_path)
polyline2 = read_ply_edge(file_path)

# Plot the polyline using matplotlib
plt.figure()
#plt.plot(polyline1[:, 0], polyline1[:, 1], marker='o', color='r', linestyle='-')
plt.plot(polyline2[:, 0], polyline2[:, 1], marker='o', color='b', linestyle='-')
plt.title("Polyline from PLY File")
plt.xlabel("X")
plt.ylabel("Y")
plt.grid(True)
plt.show()