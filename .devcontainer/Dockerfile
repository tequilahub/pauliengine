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
    ninja-build \
    && curl -L https://github.com/symengine/symengine/archive/v0.14.0.tar.gz | tar -xz \
    && rm -f v0.14.0.tar.gz \
    && cmake -S"symengine-0.14.0" -B"symengine-0.14.0/build" -DBUILD_TESTS=OFF -DBUILD_BENCHMARKS=OFF -GNinja \
    && cmake --build "symengine-0.14.0/build" --target install \
    && rm -rf symengine-0.14.0 \
    && apt-get autoremove -y \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*

COPY --from=uv --chown=vscode: /uv /uvx /bin/
