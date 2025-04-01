# Continuous Heat Method for 3D Fermat Spirals

This repository implements a C++/Python library for computing continuous heat-based Fermat spirals on 3D surfaces. The implementation combines the continuous heat method with Fermat spiral patterns to generate evenly-spaced, direction-aware spiral patterns on arbitrary 3D meshes.

![Fermat Spiral Example](data/fermatcpp.gif)

## Features

- Fast C++ implementation with Python bindings
- Support for arbitrary 3D triangle meshes in .off format
- Computation of heat-based distance fields
- Generation of Fermat spiral patterns
- Continuous heat method for smooth distance computation
- Visualization support for spiral patterns

## Installation

### Requirements

- C++ compiler with C++11 support
- CMake (>= 3.1)
- Python 3.x
- NumPy
- pybind11 (automatically installed via pip)

### Building from Source

```bash
pip install .
```

This will:
1. Build the C++ library using CMake
2. Compile the Python bindings
3. Install the package in your Python environment

## Usage

### Python API

```python
import numpy as np
from heat_fermat_3d import generate_spiral, compute_distance_field, generate_isolines

# Load your mesh vertices and faces
vertices = np.array(...)  # Nx3 array of vertex positions
faces = np.array(...)     # Mx3 array of face indices

# Generate spiral pattern
spiral_points = generate_spiral(vertices, faces, source_idx=0, num_isolines=200)

# Or compute components separately
distance_field = compute_distance_field(vertices, faces, source_idx=0)
spiral_points = generate_isolines(vertices, faces, distance_field, num_isolines=200)
```

### Example with Sample Data

```python
import numpy as np
import trimesh

# Load sample mesh
mesh = trimesh.load_mesh("data/131.off")
vertices = np.array(mesh.vertices)
faces = np.array(mesh.faces)

# Generate spiral
from heat_fermat_3d import generate_spiral
spiral_points = generate_spiral(vertices, faces)

# Visualize
mesh.vertices = vertices
mesh.show()
```

## Implementation Details

The library implements two main components:

1. **Continuous Heat Method**: Computes smooth distance fields on meshes by solving the heat equation. This provides a robust foundation for generating evenly-spaced patterns.

2. **Fermat Spiral Generation**: Uses the computed distance field to generate spiral patterns that follow the Fermat spiral equation, adapted to work on 3D surfaces.

Key C++ classes:
- `SpiralDt`: Handles distance field computation using the heat method
- `Spirals`: Manages spiral pattern generation and isoline computation

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
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.
