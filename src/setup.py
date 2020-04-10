from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

examples_extension = Extension(
	name="ctest",
	sources=["ctest.pyx"],
	libraries=["psfexlib"],
	library_dirs=[".", "fits", "levmar", "wcs", "m"],
	include_dirs=[".", "fits", "levmar", "wcs"]
)

setup(
	name="ctest",
	ext_modules = cythonize([examples_extension], compiler_directives={'language_level' : "3"})
)
