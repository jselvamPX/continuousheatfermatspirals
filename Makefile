install-ci:  ## Install the project in editable mode with test libs.
	@uv pip install -v . --system

test-import:
	@uv run python -c "import heat_fermat_3d;"