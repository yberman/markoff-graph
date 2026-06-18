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
            cmd = [
                "gcc",
                "-O3",
                "-DNDEBUG",
                "-shared",
                "-fPIC",
                "-o",
                str(out),
                str(SRC),
                "-lm",
            ]
        elif system == "Darwin":
            out = PKG / "libmarkoff.dylib"
            cmd = [
                "clang",
                "-O3",
                "-DNDEBUG",
                "-dynamiclib",
                "-o",
                str(out),
                str(SRC),
            ]
        elif system == "Windows":
            # MSVC does not have a literal /O3 flag.  /O2 is the standard
            # "maximize speed" release optimization flag for cl.exe.
            out = PKG / "libmarkoff.dll"
            cmd = [
                "cl",
                "/nologo",
                "/O2",
                "/DNDEBUG",
                "/LD",
                f"/Fe{out}",
                str(SRC),
            ]
        else:
            raise RuntimeError(f"unsupported platform: {system}")

        print("+", " ".join(cmd))
        subprocess.check_call(cmd)
        super().run()


setup(
    distclass=BinaryDistribution,
    cmdclass={"build_py": BuildPy},
)
