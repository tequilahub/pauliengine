#!/usr/bin/env bash

# install project with all dependencies
uv sync --all-groups --all-extras

# install pre-commit hooks
# uv run pre-commit install --install-hooks
