install-ci:  ## Install the project in editable mode with test libs.
	@uv venv
	@uv run python -m pip install -v .

test-import:
	@uv run python -c "import heat_fermat_3d;"