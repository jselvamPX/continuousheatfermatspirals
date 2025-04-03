import os
import sys
import importlib.util

# Try to import from site-packages first
try:
    from . import _heat_fermat_3d_core
except ImportError:
    # If that fails, try to find the module in site-packages
    for site_dir in sys.path:
        module_path = os.path.join(
            site_dir, "heat_fermat_3d", "_heat_fermat_3d_core.so"
        )
        if os.path.exists(module_path):
            spec = importlib.util.spec_from_file_location(
                "_heat_fermat_3d_core", module_path
            )
            _heat_fermat_3d_core = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(_heat_fermat_3d_core)
            break
    else:
        raise ImportError(
            "Could not find the _heat_fermat_3d_core extension module.\n"
            f"Current Python interpreter: {sys.executable}\n"
            f"Current working directory: {os.getcwd()}\n"
            f"Site-packages directories: {[p for p in sys.path if 'site-packages' in p]}\n"
            f"Package directory: {os.path.dirname(__file__)}\n"
            f"Files in package directory: {os.listdir(os.path.dirname(__file__))}\n"
            "Please make sure the package is installed correctly in the current Python environment."
        )

# Import the functions from the extension module
generate_spiral = _heat_fermat_3d_core.generate_spiral
compute_distance_field = _heat_fermat_3d_core.compute_distance_field
generate_isolines = _heat_fermat_3d_core.generate_isolines

__version__ = "0.1.0"
__all__ = ["generate_spiral", "compute_distance_field", "generate_isolines"]
