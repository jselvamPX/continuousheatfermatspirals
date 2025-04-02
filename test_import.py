import sys
import os
import importlib

# Print Python information
print(f"Python executable: {sys.executable}")
print(f"Python version: {sys.version}")
print(f"Current working directory: {os.getcwd()}")

# Try to import the package
try:
    import heat_fermat_3d

    print("Import successful!")
    print(f"Package location: {heat_fermat_3d.__file__}")
    print(f"Available functions: {dir(heat_fermat_3d)}")
except ImportError as e:
    print(f"Import failed: {e}")

# Try to find the module directly
try:
    import site

    for site_dir in site.getsitepackages():
        package_dir = os.path.join(site_dir, "heat_fermat_3d")
        if os.path.exists(package_dir):
            print(f"Found package directory: {package_dir}")
            print(f"Files in package directory: {os.listdir(package_dir)}")
except ImportError:
    pass
