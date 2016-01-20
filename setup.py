# setup.py --- setup routine for wrapping C++ H corona code

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

setup(ext_modules = cythonize(Extension(
    "py_los_profile",
    ["py_los_profile.pyx"],                 # our Cython source
    language="c++",             # generate C++ code
    extra_compile_args=["-O3"],
    extra_link_args=["-lm"],
    define_macros=[('SRCFNSLOC','"./"')]))
)
