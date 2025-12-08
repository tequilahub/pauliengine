# PauliEngine: Fast Arithmetics for Quantum Operators

This project is currently in **beta**. Some features are incomplete or still under development.  
The following items are planned or not yet fully implemented:

-   A detailed online Documentation of PauliEngines functionality.
    
-   Installation via pip.
    
-   A comprehensive test suite.

# Example I 

### Operator Construction
```python
import pauliengine as pe
# Example : Z on qubit 0 with coefficient 1.0
p1 = pe.PauliString ((1.0, {0: "Z"}))

# Example : X on qubit 1 with coefficient "a" ( symbolic )
p2 = pe.PauliString (("a", {1: "X"}))

# Example : OpenFermion style
p3 = pe.PauliString ((1.0, [("X", 0), ("Y", 2) ]))
```
### Standard Operations
```python
# Addition:
p4 = p2 + p3
# Multiplication
p5 = p4 * p1
# fast commutators
c = p1.commutator(p2)
```

### Parametrized Operators
```python
# create PauliString objects
p1 = pe.PauliString ((1.0 , {0: "Z"})) # not parametrized
p2 = pe.PauliString (("a", {1: "X"}))  # parametrized
# assemble operator
H = QubitHamiltonian([p1 , p2 ])
```
In beta version, mixed types need to be created like this (create atomic PauliStrings, then assemble operator)

### Operations on Parametrized Operators
```python
# Differentiate
dH = H.diff ("a")
# Substitute
H2 = H.subs ({"a": 2.0})
```

# Installation and Build Guide 

  

This guide explains how to install and build **SymEngine** and integrate it into a C++ project on **Windows** using **Conan**, **CMake**, and **Visual Studio**.

  

---

  

## Prerequisites

  

Before you begin, make sure you have the following installed:

  

-  **Git** (installed and accessible from the command line)

  

-  **CMake** (installed and added to your PATH)

  

-  **Python** (for installing Conan)

 


## Step 1: Install Conan

  

Install  Conan  via  `pip`:
```bash
pip  install  conan

conan  profile  detect

```

---

## Step 2: Install SymEngine

  MacOS only
  ```bash
conan install --build=symengine/0.14.0
```
## Step 3 Clone Project
Clone the project into desired directory.
If nanobind is missing, clone nanobind with the option `--recursive` into `src`.

---

## Step 4: Integrate SymEngine in Your Project
In the directory `src`
  
```bash
conan  install  .  --output-folder=build  --build=missing
```

---

## Step 5: Configure with CMake

  ```bash

cd  build

cmake  ..  -DCMAKE_TOOLCHAIN_FILE=conan_toolchain.cmake
```

---

## Step 6: Build the Project

  ```bash

cmake  --build  .  --config  Release
```
---
## Step 7: Using the module
When building this project as a Python extension module using **nanobind**, the resulting file type depends on the operating system:

-   **Linux:** `.so`
    
-   **macOS:** `.so`
    
-   **Windows:** `.pyd`
    

Although Linux and macOS both generate `.so` files, they use different binary formats (ELF on Linux, Mach-O on macOS). On Windows, Python loads native modules using the `.pyd` extension, which is similar to a DLL.

The final compiled extension module will be located in:
```bash
src/build/Release
```
This file can now be used in any Python project, and its functionality can be accessed via:

```bash
 import PauliEngine
```




