from __future__ import annotations

import os
import platform
import subprocess
from pathlib import Path

from setuptools import setup
from setuptools.command.build_py import build_py


ROOT = Path(__file__).parent
PKG = ROOT / "src" / "markoff_graph"
SRC = PKG / "libmarkoff.c"


class BuildPy(build_py):
    def run(self):
        self.build_native_library()
        super().run()

    def build_native_library(self):
        system = platform.system()

        if system == "Linux":
            out = PKG / "libmarkoff.so"
            cmd = ["gcc", "-O3", "-shared", "-fPIC", "-o", str(out), str(SRC)]

        elif system == "Darwin":
            out = PKG / "libmarkoff.dylib"
            cmd = ["clang", "-O3", "-dynamiclib", "-o", str(out), str(SRC)]

        elif system == "Windows":
            out = PKG / "libmarkoff.dll"
            cc = os.environ.get("CC", "gcc")
            cmd = [cc, "-O3", "-shared", "-o", str(out), str(SRC)]

        else:
            raise RuntimeError(f"Unsupported platform: {system}")

        print("+", " ".join(cmd))
        subprocess.check_call(cmd)


setup(cmdclass={"build_py": BuildPy})
