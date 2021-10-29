"""
Setup file to Cythonize FASTQ Reader using "python3 setup.py build_ext --inplace"
"""
import os
from distutils.core import setup
from Cython.Build import cythonize

# reader_file = '{0}{1}FileWriter.pyx'.format(os.path.dirname(__file__), os.sep)
# reader_file = '{0}{1}FASTQReader.pyx'.format(os.path.dirname(__file__), os.sep)

setup(
    name="Valkyries",
    author='Dennis Simpson',
    author_email='dennis@email.unc.edu',
    ext_modules=cythonize(["FileWriter.pyx", "FASTQReader.pyx"], annotate=False)
)

