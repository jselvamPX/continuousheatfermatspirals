#!/bin/bash

# This script helps diagnose and fix issues with the heat_fermat_3d package

# Print the current Python environment
echo "Current Python interpreter:"
which python
echo "Python version:"
python --version

# Uninstall any existing installations
echo -e "\nUninstalling any existing heat_fermat_3d package..."
pip uninstall -y heat_fermat_3d

# Install the package
echo -e "\nInstalling heat_fermat_3d package..."
pip install -v .

# Check the installation
echo -e "\nChecking installation..."
pip show heat_fermat_3d

# List the installed files
echo -e "\nListing installed files:"
pip_install_dir=$(pip show heat_fermat_3d | grep Location | awk '{print $2}')
echo "Package installed in: $pip_install_dir"
echo "Files in heat_fermat_3d directory:"
ls -la "$pip_install_dir/heat_fermat_3d/"

# Test importing the package
echo -e "\nTesting import..."
python -c "
import os
import sys
print(f'Python executable: {sys.executable}')
print(f'Python version: {sys.version}')
print(f'Python path: {sys.path}')
print('Attempting to import heat_fermat_3d...')
try:
    import heat_fermat_3d
    print('Import successful!')
    print(f'Available functions: {dir(heat_fermat_3d)}')
except ImportError as e:
    print(f'Import failed: {e}')
"

echo -e "\nSetup and test complete."
