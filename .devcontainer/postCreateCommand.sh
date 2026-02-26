#!/usr/bin/env bash

set -euxo pipefail

# install uv and uvx shell completions
echo 'eval "$(uv generate-shell-completion zsh)"' >> /home/vscode/.zshrc
echo 'eval "$(uvx --generate-shell-completion zsh)"' >> /home/vscode/.zshrc

# set up GDB to view contents of STL containers
gcc_version="$(gcc --version | awk 'NR==1 {print $NF}')"
cat <<EOT >> ~/.gdbinit
python
import sys
sys.path.insert(0, "/usr/share/gcc-$gcc_version/python")
from libstdcxx.v6.printers import register_libstdcxx_printers
register_libstdcxx_printers (None)
end
EOT

# install conan
uv venv
source /home/vscode/.venv/bin/activate
uv pip install conan 

# setup conan profile and install dependencies
conan profile detect --force 
conan install .  --output-folder=build --build=missing -pr ./conan_profile