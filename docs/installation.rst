Installation
============

From PyPI
---------

.. code-block:: bash

   pip install pauliengine

Prebuilt wheels are available for Linux (x86_64/aarch64), macOS (arm64) and
Windows (x86_64/arm64) on Python 3.10–3.13.

From source
-----------

Building from source needs a C++20 compiler (MSVC 19.3+, GCC 11+ or Clang 14+)
and CMake 3.20+. There are no C++ libraries to install beforehand:

.. code-block:: bash

   pip install .

Use ``pip install -e .`` for an editable development install.

Building this documentation
---------------------------

With PauliEngine installed, from the repository root:

.. code-block:: bash

   pip install sphinx breathe furo
   sphinx-build -b html docs docs/_build/html

The C++ reference additionally needs `Doxygen <https://www.doxygen.nl>`_ on
``PATH``; without it, only the Python reference is built.
