#!/usr/bin/env bash

set -euo pipefail

# Script to install dependencies for pauliengine project
# Usage: ./install-deps.sh [install_prefix] [--skip-boost-test] [--skip-msgpack] [--help]

show_help() {
    cat << EOF
Usage: $0 [INSTALL_PREFIX] [OPTIONS]

Install C++ dependencies for pauliengine project.

Arguments:
    INSTALL_PREFIX      Directory to install dependencies (default: /usr/local)

Options:
    --help, -h          Show this help message

Examples:
    $0                                   # Install all deps to /usr/local
    $0 \$HOME/Software                   # Install all deps to \$HOME/Software
EOF
}

# Default values
DEFAULT_PREFIX="/usr/local"
INSTALL_PREFIX="$DEFAULT_PREFIX"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --help|-h)
            show_help
            exit 0
            ;;
        --*)
            echo "Unknown option: $1" >&2
            echo "Use --help for usage information." >&2
            exit 1
            ;;
        *)
            # Assume it's the install prefix if it doesn't start with --
            if [[ -z "${INSTALL_PREFIX_SET:-}" ]]; then
                INSTALL_PREFIX="$1"
                INSTALL_PREFIX_SET=true
            else
                echo "Multiple install prefixes specified: $INSTALL_PREFIX and $1" >&2
                exit 1
            fi
            shift
            ;;
    esac
done

echo "Installing C++ dependencies to: $INSTALL_PREFIX"

# Create install directory if it doesn't exist
mkdir -p "$INSTALL_PREFIX"


curl -L https://github.com/symengine/symengine/archive/master.tar.gz | tar -xz
rm -f master.tar.gz

cmake -S"symengine-master" -B"symengine-master/build" \
  -DCMAKE_INSTALL_PREFIX="$INSTALL_PREFIX" \
  -DBUILD_TESTS=OFF \
  -DBUILD_BENCHMARKS=OFF \
  -GNinja

cmake --build "symengine-master/build" --target install
rm -rf symengine-master

if [[ "$INSTALL_PREFIX" == "$DEFAULT_PREFIX" ]]; then
    echo "Remember to export SymEngine_DIR=$INSTALL_PREFIX/lib/cmake/SymEngine"
fi

echo
echo "Dependencies installation completed successfully!"
echo "Install location: $INSTALL_PREFIX"
echo "Make sure to set CMAKE_PREFIX_PATH=$INSTALL_PREFIX when building pauliengine"
