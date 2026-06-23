"""Setup script for SiteGA Python bindings (pybind11).

Build the extension with:
    pip install .
or:
    python setup.py build_ext --inplace

The resulting module is importable as `import sitega` and provides:
    - sitega.TrainParams: parameter struct for training
    - sitega.train(params): train a model, returns (rc, mat_path, loc_path, best_fit)
"""
import os
from setuptools import setup

# Try pybind11; fall back to raw extension if not installed
try:
    from pybind11.setup_helpers import Pybind11Extension, build_ext
    ext_modules = [
        Pybind11Extension(
            "sitega",
            sources=["src/bindings.cpp", "src/andy05cell.cpp"],
            include_dirs=["src"],
            extra_compile_args=[
                "-O3", "-ffast-math", "-funroll-loops",
                "-fopenmp",
            ],
            extra_link_args=["-fopenmp"],
            cxx_std=17,
        )
    ]
    cmdclass = {"build_ext": build_ext}
except ImportError:
    from setuptools import Extension
    ext_modules = [
        Extension(
            "sitega",
            sources=["src/bindings.cpp", "src/andy05cell.cpp"],
            include_dirs=["src"],
            extra_compile_args=[
                "-O3", "-ffast-math", "-funroll-loops",
                "-fopenmp",
            ],
            extra_link_args=["-fopenmp"],
        )
    ]
    cmdclass = {}

setup(
    name="sitega",
    version="1.0.0",
    ext_modules=ext_modules,
    cmdclass=cmdclass,
    zip_safe=False,
)