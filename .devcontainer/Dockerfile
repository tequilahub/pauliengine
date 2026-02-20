FROM ghcr.io/astral-sh/uv:latest AS uv

FROM mcr.microsoft.com/devcontainers/base:2-ubuntu24.04

# install needed packages
RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
    && apt-get -y install \
    build-essential \
    cmake \
    gcc \
    gdb \
    libgmp-dev \
    ninja-build


COPY --from=uv --chown=vscode /uv /uvx /bin/
