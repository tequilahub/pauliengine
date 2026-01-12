FROM ghcr.io/astral-sh/uv:latest AS uv

FROM mcr.microsoft.com/devcontainers/base:2-ubuntu24.04

# Install needed packages, and only install MKL on x86_64
RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
    && apt-get -y install \
    build-essential \
    cmake \
    gcc \
    gdb \
    libblis-openmp-dev \
    libopenblas-openmp-dev \
    ninja-build \
    openmpi-bin \
    && if [ "$(uname -m)" = "x86_64" ]; then \
         apt-get -y install wget gpg-agent && \
         wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | gpg --dearmor | tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null && \
         echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" > /etc/apt/sources.list.d/oneAPI.list && \
         apt-get update && \
         apt-get -y install intel-oneapi-mkl-core-devel && \
         apt-get -y remove gpg-agent wget; \
       fi \
    && apt-get autoremove -y \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*

COPY --from=uv --chown=vscode: /uv /uvx /bin/