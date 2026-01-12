#!/usr/bin/env bash

# source OneAPI environment
. /opt/intel/oneapi/setvars.sh &> /dev/null

# install project with all dependencies
# uv sync --all-groups --all-extras

# install pre-commit hooks
# uv run pre-commit install --install-hooks