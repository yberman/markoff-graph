from pathlib import Path
import platform
import subprocess

from setuptools import Distribution, setup
from setuptools.command.build_py import build_py


ROOT = Path(__file__).parent
PKG = ROOT / "src" / "markoff_graph"
SRC = PKG / "libmarkoff.c"


class BinaryDistribution(Distribution):
    def has_ext_modules(self):
        return True


class BuildPy(build_py):
    def run(self):
        system = platform.system()

        if system == "Linux":
            out = PKG / "libmarkoff.so"
            cmd = ["gcc", "-O3", "-shared", "-fPIC", "-o", str(out), str(SRC)]
        elif system == "Darwin":
            out = PKG / "libmarkoff.dylib"
            cmd = ["clang", "-O3", "-dynamiclib", "-o", str(out), str(SRC)]
        elif system == "Windows":
            out = PKG / "libmarkoff.dll"
            cmd = ["cl", "/nologo", "/O2", "/LD", f"/Fe{out}", str(SRC)]
        else:
            raise RuntimeError(f"unsupported platform: {system}")

        print("+", " ".join(cmd))
        subprocess.check_call(cmd)
        super().run()


setup(
    distclass=BinaryDistribution,
    cmdclass={"build_py": BuildPy},
)
