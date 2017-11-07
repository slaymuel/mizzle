from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

ext_module = Extension(
    "potential",
    ["drizzle/potential.pyx"],
    extra_compile_args = ["-O3", "-ffast-math", "-march=native", "-fopenmp" ],
    extra_link_args=['-fopenmp'],
    include_dirs=[numpy.get_include()],
)

setup(name='drizzle',
      version='1.0',
      description='Wet arbitrary metal-oxide surfaces',
      author='Samuel Stenberg',
      author_email='samuel.stenberg@hotmail.com',
      license='MIT',
      packages=['drizzle'],
      install_requires=['argcomplete','numpy','pandas','mdtraj','tqdm', 'scipy', 'radish'],
      include_package_data=True,
      cmdclass = {'build_ext': build_ext},
      ext_modules = [ext_module],
      zip_safe=False)