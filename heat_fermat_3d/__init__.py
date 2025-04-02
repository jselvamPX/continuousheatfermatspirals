import os
import sys
import importlib.util
from pathlib import Path

# Get the directory where this __init__.py file is located
package_dir = Path(__file__).parent.absolute()

# Try to find the extension module
extension_path = None
for file in os.listdir(package_dir):
    if file.startswith("_heat_fermat_3d_core") and file.endswith(".so"):
        extension_path = os.path.join(package_dir, file)
        break

if extension_path is None:
    # If we're running from the source directory, try to find the module in the site-packages
    try:
        import site

        for site_dir in site.getsitepackages():
            package_site_dir = os.path.join(site_dir, "heat_fermat_3d")
            if os.path.exists(package_site_dir):
                for file in os.listdir(package_site_dir):
                    if file.startswith("_heat_fermat_3d_core") and file.endswith(".so"):
                        extension_path = os.path.join(package_site_dir, file)
                        break
                if extension_path is not None:
                    break
    except ImportError:
        pass

if extension_path is None:
    raise ImportError(
        "Could not find the _heat_fermat_3d_core extension module. "
        "Please make sure the package is installed correctly."
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
