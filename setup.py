from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy
from os import path
from glob import glob


setup(
    name="eiliftover",
    ext_modules = cythonize([Extension(path.join("eiliftover.util.contrast"),
                                     [path.join("eiliftover", "util", "contrast.pyx")]),
                             Extension(path.join("eiliftover.util.overlap"),
                                       [path.join("eiliftover", "util", "overlap.pyx")])
                             ]),
    include_dirs=[numpy.get_include()],
    scripts=glob("util/*.py")
)