import numpy as np
import trimesh
from heat_fermat_3d import generate_spiral, compute_distance_field, generate_isolines

# Load sample mesh
mesh = trimesh.load_mesh("data/131.off")
vertices = np.array(mesh.vertices)
faces = np.array(mesh.faces)

# Generate spiral
print("Generating spiral...")
spiral_points = generate_spiral(vertices, faces)
print(f"Generated spiral with {len(spiral_points)} points")

# Or compute components separately
print("\nComputing distance field...")
distance_field = compute_distance_field(vertices, faces, source_idx=0)
print(f"Computed distance field with {len(distance_field)} values")

print("\nGenerating isolines...")
spiral_points = generate_isolines(vertices, faces, distance_field, num_isolines=200)
print(f"Generated isolines with {len(spiral_points)} points")

# Save the spiral points to a file
np.savetxt("spiral_points.txt", spiral_points)
print("\nSaved spiral points to spiral_points.txt")
