# Fermat3D Python Bindings

Python bindings for the Fermat3D library, which generates continuous spiral lines on 3D meshes using the heat method.

## Overview

The Python bindings expose the following functions:

- `generate_spiral(vertices, faces, source_idx=0, num_isolines=200)`: Generate a continuous spiral line from vertices and faces
- `compute_distance_field(vertices, faces, source_idx=0)`: Compute a distance field from vertices and faces
- `generate_isolines(vertices, faces, distance_field, num_isolines=200)`: Generate isolines from vertices, faces, and a distance field

## Current Implementation Status

The Python bindings are now fully functional and use the C++ implementation to generate the spiral. The bindings are implemented using pybind11 and provide a high-level API to the Fermat3D library.

## Usage

Here's a simple example of how to use the Python module:

```python
import numpy as np
import os
import sys
import importlib.util

# Import the C++ extension directly
extension_path = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "heat_fermat_3d.cpython-311-x86_64-linux-gnu.so",
)

# Import the C++ extension
spec = importlib.util.spec_from_file_location("heat_fermat_3d", extension_path)
heat_fermat_3d = importlib.util.module_from_spec(spec)
spec.loader.exec_module(heat_fermat_3d)

# Load vertices and faces (e.g., from an OBJ or OFF file)
vertices = np.array([...])  # Nx3 array of vertex positions
faces = np.array([...])     # Mx3 array of face indices

# Generate the spiral
spiral_points = heat_fermat_3d.generate_spiral(vertices, faces, source_idx=0, num_isolines=200)

# Save or visualize the spiral
# ...
```

For a complete example, see the `example.py` file in the repository.

## Example

The `example.py` file demonstrates how to:

1. Load a mesh from an OFF file
2. Generate a spiral using the `heat_fermat_3d` module
3. Save the spiral to an OBJ file
4. Visualize the mesh and spiral (if trimesh and matplotlib are available)

To run the example:

```bash
# Run the example with the default mesh
python example.py

# Or specify a custom mesh file
python example.py path/to/your/mesh.off
```

## API Reference

### generate_spiral

```python
def generate_spiral(vertices, faces, source_idx=0, num_isolines=200):
    """
    Generate a continuous spiral line from vertices and faces.
    
    Parameters:
    -----------
    vertices : numpy.ndarray
        Nx3 array of vertex positions
    faces : numpy.ndarray
        Mx3 array of face indices
    source_idx : int, optional
        Index of the source vertex for the heat method (default: 0)
    num_isolines : int, optional
        Number of isolines to generate (default: 200)
    
    Returns:
    --------
    numpy.ndarray
        Px3 array of points along the spiral
    """
```

### compute_distance_field

```python
def compute_distance_field(vertices, faces, source_idx=0):
    """
    Compute a distance field from vertices and faces using the heat method.
    
    Parameters:
    -----------
    vertices : numpy.ndarray
        Nx3 array of vertex positions
    faces : numpy.ndarray
        Mx3 array of face indices
    source_idx : int, optional
        Index of the source vertex for the heat method (default: 0)
    
    Returns:
    --------
    numpy.ndarray
        N-length array of distance values for each vertex
    """
```

### generate_isolines

```python
def generate_isolines(vertices, faces, distance_field, num_isolines=200):
    """
    Generate isolines from vertices, faces, and a distance field.
    
    Parameters:
    -----------
    vertices : numpy.ndarray
        Nx3 array of vertex positions
    faces : numpy.ndarray
        Mx3 array of face indices
    distance_field : numpy.ndarray
        N-length array of distance values for each vertex
    num_isolines : int, optional
        Number of isolines to generate (default: 200)
    
    Returns:
    --------
    numpy.ndarray
        Px3 array of points along the spiral
    """
