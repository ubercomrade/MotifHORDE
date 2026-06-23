from __future__ import annotations

import os
import sys
from distutils.errors import CCompilerError, CompileError, LinkError
from pathlib import Path

from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

SITEGA_DIR = Path("src") / "sitega"


def _compiler_command(compiler) -> str:
    for name in ("compiler_cxx", "compiler_so", "compiler"):
        command = getattr(compiler, name, None)
        if command:
            if isinstance(command, (list, tuple)):
                return " ".join(command)
            return str(command)
    return ""


def _is_gcc(command: str) -> bool:
    if not command:
        return False
    compiler_name = os.path.basename(command.split()[0]).lower()
    command = command.lower()
    return (
        "gcc" in compiler_name or "g++" in compiler_name
    ) and "clang" not in command


def _use_openmp(compiler_type: str, command: str) -> bool:
    return (
        sys.platform.startswith("linux")
        and compiler_type == "unix"
        and _is_gcc(command)
    )


def _compile_args(compiler_type: str, openmp: bool) -> list[str]:
    if compiler_type == "msvc":
        return ["/O2"]

    args = ["-O3"]
    if openmp:
        args.append("-fopenmp")
    return args


def _link_args(compiler_type: str, openmp: bool) -> list[str]:
    if compiler_type == "msvc":
        return []
    if openmp:
        return ["-fopenmp"]
    return []


class BuildExt(build_ext):
    def build_extension(self, ext) -> None:
        compiler_type = self.compiler.compiler_type
        command = _compiler_command(self.compiler)
        openmp = _use_openmp(compiler_type, command)
        self._set_build_args(ext, compiler_type, openmp)

        if not openmp:
            super().build_extension(ext)
            return

        try:
            super().build_extension(ext)
        except (CCompilerError, CompileError, LinkError):
            self.announce(
                "OpenMP build failed; retrying SiteGA without OpenMP",
                level=3,
            )
            self._set_build_args(ext, compiler_type, openmp=False)
            old_force = self.force
            self.force = True
            try:
                super().build_extension(ext)
            finally:
                self.force = old_force

    @staticmethod
    def _set_build_args(ext, compiler_type: str, openmp: bool) -> None:
        ext.extra_compile_args = _compile_args(compiler_type, openmp)
        ext.extra_link_args = _link_args(compiler_type, openmp)


ext_modules = [
    Pybind11Extension(
        "sitega",
        sources=[
            str(SITEGA_DIR / "bindings.cpp"),
            str(SITEGA_DIR / "andy05cell.cpp"),
        ],
        include_dirs=[str(SITEGA_DIR)],
        cxx_std=17,
    )
]


setup(
    ext_modules=ext_modules,
    cmdclass={"build_ext": BuildExt},
    zip_safe=False,
)
