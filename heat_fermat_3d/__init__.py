import os
import sys
import importlib.util
from pathlib import Path
import glob

# Get the directory where this __init__.py file is located
package_dir = Path(__file__).parent.absolute()

# Try to find the extension module using multiple approaches
extension_path = None

# Approach 1: Look for the module in the same directory as this file
for file in os.listdir(package_dir):
    if file.startswith("_heat_fermat_3d_core") and (
        file.endswith(".so") or file.endswith(".pyd")
    ):
        extension_path = os.path.join(package_dir, file)
        break

# Approach 2: Use glob to find the module with any extension
if extension_path is None:
    possible_paths = glob.glob(os.path.join(package_dir, "_heat_fermat_3d_core*"))
    if possible_paths:
        extension_path = possible_paths[0]

# Approach 3: Look in site-packages
if extension_path is None:
    try:
        import site

        for site_dir in site.getsitepackages():
            package_site_dir = os.path.join(site_dir, "heat_fermat_3d")
            if os.path.exists(package_site_dir):
                # Try direct filename
                for file in os.listdir(package_site_dir):
                    if file.startswith("_heat_fermat_3d_core") and (
                        file.endswith(".so") or file.endswith(".pyd")
                    ):
                        extension_path = os.path.join(package_site_dir, file)
                        break

                # Try glob
                if extension_path is None:
                    possible_paths = glob.glob(
                        os.path.join(package_site_dir, "_heat_fermat_3d_core*")
                    )
                    if possible_paths:
                        extension_path = possible_paths[0]

                if extension_path is not None:
                    break
    except (ImportError, AttributeError):
        pass

# Approach 4: Try to directly import using ctypes as a last resort
if extension_path is None:
    try:
        import ctypes

        try:
            # Try Linux/macOS naming
            _heat_fermat_3d_core = ctypes.CDLL(
                os.path.join(package_dir, "_heat_fermat_3d_core.so")
            )
            extension_path = os.path.join(package_dir, "_heat_fermat_3d_core.so")
        except OSError:
            try:
                # Try Windows naming
                _heat_fermat_3d_core = ctypes.CDLL(
                    os.path.join(package_dir, "_heat_fermat_3d_core.pyd")
                )
                extension_path = os.path.join(package_dir, "_heat_fermat_3d_core.pyd")
            except OSError:
                pass
    except ImportError:
        pass

# If still not found, raise a detailed error
if extension_path is None:
    # Provide more detailed error information
    python_path = sys.executable
    current_dir = os.getcwd()
    site_packages = []
    try:
        import site

        site_packages = site.getsitepackages()
    except ImportError:
        pass

    raise ImportError(
        f"Could not find the _heat_fermat_3d_core extension module.\n"
        f"Current Python interpreter: {python_path}\n"
        f"Current working directory: {current_dir}\n"
        f"Site-packages directories: {site_packages}\n"
        f"Package directory: {package_dir}\n"
        f"Files in package directory: {os.listdir(package_dir)}\n"
        "Please make sure the package is installed correctly in the current Python environment."
    )

# Load the extension module
spec = importlib.util.spec_from_file_location("_heat_fermat_3d_core", extension_path)
_heat_fermat_3d_core = importlib.util.module_from_spec(spec)
sys.modules["heat_fermat_3d._heat_fermat_3d_core"] = _heat_fermat_3d_core
spec.loader.exec_module(_heat_fermat_3d_core)

# Import the functions from the extension module
generate_spiral = _heat_fermat_3d_core.generate_spiral
compute_distance_field = _heat_fermat_3d_core.compute_distance_field
generate_isolines = _heat_fermat_3d_core.generate_isolines

__version__ = "0.1.0"
__all__ = ["generate_spiral", "compute_distance_field", "generate_isolines"]
