"""Python bindings for Fermat3D spiral generation"""

import os
import sys
import numpy as np

import os
import sys
import importlib.util


def _find_module():
    # First try local directory
    package_dir = os.path.dirname(__file__)
    for file in os.listdir(package_dir):
        if file.startswith("_heat_fermat_3d_core") and file.endswith(
            (".so", ".dylib", ".dll")
        ):
            return os.path.join(package_dir, file)

    # Then try site-packages
    try:
        import site

        for site_dir in site.getsitepackages():
            ext_dir = os.path.join(site_dir, "heat_fermat_3d")
            if os.path.exists(ext_dir):
                for file in os.listdir(ext_dir):
                    if file.startswith("_heat_fermat_3d_core") and file.endswith(
                        (".so", ".dylib", ".dll")
                    ):
                        return os.path.join(ext_dir, file)
    except Exception:
        pass

    return None


_module_path = _find_module()
if _module_path is None:
    raise ImportError(
        "Could not find heat_fermat_3d extension module. "
        "This could be because the package was not built correctly. "
        f"Searched in {os.path.dirname(__file__)} and site-packages"
    )

try:
    _spec = importlib.util.spec_from_file_location("_heat_fermat_3d_core", _module_path)
    _core = importlib.util.module_from_spec(_spec)
    sys.modules["_heat_fermat_3d_core"] = _core
    _spec.loader.exec_module(_core)
except Exception as e:
    raise ImportError(
        f"Failed to load heat_fermat_3d extension module from {_module_path}. "
        f"Error: {e}"
    ) from e


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
    vertices = np.asarray(vertices, dtype=np.float64)
    faces = np.asarray(faces, dtype=np.int32)

    if vertices.ndim != 2 or vertices.shape[1] != 3:
        raise ValueError("vertices must be an Nx3 array")
    if faces.ndim != 2 or faces.shape[1] != 3:
        raise ValueError("faces must be an Mx3 array")
    if source_idx < 0 or source_idx >= len(vertices):
        raise ValueError("source_idx must be a valid vertex index")
    if num_isolines < 1:
        raise ValueError("num_isolines must be positive")

    return _core.generate_spiral(vertices, faces, source_idx, num_isolines)


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
    vertices = np.asarray(vertices, dtype=np.float64)
    faces = np.asarray(faces, dtype=np.int32)

    if vertices.ndim != 2 or vertices.shape[1] != 3:
        raise ValueError("vertices must be an Nx3 array")
    if faces.ndim != 2 or faces.shape[1] != 3:
        raise ValueError("faces must be an Mx3 array")
    if source_idx < 0 or source_idx >= len(vertices):
        raise ValueError("source_idx must be a valid vertex index")

    return _core.compute_distance_field(vertices, faces, source_idx)


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
    vertices = np.asarray(vertices, dtype=np.float64)
    faces = np.asarray(faces, dtype=np.int32)
    distance_field = np.asarray(distance_field, dtype=np.float64)

    if vertices.ndim != 2 or vertices.shape[1] != 3:
        raise ValueError("vertices must be an Nx3 array")
    if faces.ndim != 2 or faces.shape[1] != 3:
        raise ValueError("faces must be an Mx3 array")
    if distance_field.ndim != 1 or len(distance_field) != len(vertices):
        raise ValueError(
            "distance_field must be a 1D array with length equal to number of vertices"
        )
    if num_isolines < 1:
        raise ValueError("num_isolines must be positive")

    return _core.generate_isolines(vertices, faces, distance_field, num_isolines)


# Define the public API
__version__ = "0.1.0"
__all__ = ["generate_spiral", "compute_distance_field", "generate_isolines"]
