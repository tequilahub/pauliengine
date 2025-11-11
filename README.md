# ⚙️ SymEngine Installation and Build Guide (Windows)

This guide explains how to install and build **SymEngine** and integrate it into a C++ project on **Windows** using **Conan**, **CMake**, and **Visual Studio**.

---

## 🧰 Prerequisites

Before you begin, make sure you have the following installed:

- **Visual Studio** (with C++ build tools)  
  During installation, select:
  - ✅ *Desktop development with C++*

- **Git** (installed and accessible from the command line)

- **CMake** (installed and added to your PATH)

- **Python** (for installing Conan)

---
```bash
##  Step 1: Install Conan

Install Conan via `pip`:


pip install conan
---
##  Step 2: Install SymEngine

conan profile detect
---
## Step 3: Integrate SymEngine in Your Project

conan install . --output-folder=build --build=missing
---
## Step 4: Configure with CMake

cd build
cmake .. -DCMAKE_TOOLCHAIN_FILE=conan_toolchain.cmake
---
## Step 5: Build the Project

cmake --build . --config Release
---