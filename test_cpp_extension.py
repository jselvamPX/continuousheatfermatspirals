import numpy as np
import heat_fermat_3d


def load_off_file(file_path):
    """Load vertices and faces from an OFF file."""
    with open(file_path, "r") as f:
        lines = f.readlines()

    # Skip comments
    line_idx = 0
    while lines[line_idx].startswith("#"):
        line_idx += 1

    # Check if the file is an OFF file
    if not lines[line_idx].strip() == "OFF":
        raise ValueError("File is not in OFF format")
    line_idx += 1

    # Read the number of vertices, faces, and edges
    counts = lines[line_idx].strip().split()
    num_vertices = int(counts[0])
    num_faces = int(counts[1])
    line_idx += 1

    # Read vertices
    vertices = np.zeros((num_vertices, 3))
    for i in range(num_vertices):
        vertices[i] = np.array([float(x) for x in lines[line_idx].strip().split()])
        line_idx += 1

    # Read faces
    faces = np.zeros((num_faces, 3), dtype=np.int32)
    for i in range(num_faces):
        face_data = [int(x) for x in lines[line_idx].strip().split()]
        if face_data[0] != 3:
            raise ValueError("Only triangular faces are supported")
        faces[i] = face_data[1:4]
        line_idx += 1

    return vertices, faces


# Save the spiral to an OBJ file
def save_obj_file(file_path, vertices):
    """Save vertices as a line in an OBJ file."""
    with open(file_path, "w") as f:
        for i, v in enumerate(vertices):
            f.write(f"v {v[0]} {v[1]} {v[2]}\n")

        # Add lines connecting consecutive vertices
        for i in range(len(vertices) - 1):
            f.write(f"l {i + 1} {i + 2}\n")


# Load the mesh
file_path = "data/131.off"
print(f"Loading mesh from {file_path}...")
vertices, faces = load_off_file(file_path)
print(f"Loaded mesh with {len(vertices)} vertices and {len(faces)} faces.")

# Generate the spiral
print("Generating spiral...")
spiral_points = heat_fermat_3d.generate_spiral(
    vertices, faces, source_idx=0, num_isolines=200
)
print(f"Generated spiral with {len(spiral_points)} points.")

# Save the spiral to an OBJ file
output_file = "spiral_cpp.obj"
save_obj_file(output_file, spiral_points)
print(f"Saved spiral to {output_file}")

print("Done!")
