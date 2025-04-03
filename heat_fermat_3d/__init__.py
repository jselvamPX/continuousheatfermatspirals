from . import _heat_fermat_3d_core


# Import the functions from the extension module
generate_spiral = _heat_fermat_3d_core.generate_spiral
compute_distance_field = _heat_fermat_3d_core.compute_distance_field
generate_isolines = _heat_fermat_3d_core.generate_isolines

__version__ = "0.1.0"
__all__ = ["generate_spiral", "compute_distance_field", "generate_isolines"]
