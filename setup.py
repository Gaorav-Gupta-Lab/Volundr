
from setuptools import setup, find_namespace_packages
import sys

if sys.version_info < (3, 5):
    sys.stdout.write("At least Python 3.5 is required.\n")
    sys.exit(1)
# packages=['volundr', 'Valkyries'],

setup(
    name='Volundr',
    version='1.0.0',
    packages=find_namespace_packages(include=['Volundr.*']),

    url='',
    license='MIT',
    author='Dennis Simpson',
    author_email='dennis@email.unc.edu',
    long_description=open('README.md').read(),

    description='Package for using CRISPR to do synthetic lethal and viable assays in mammalian cells',
    install_requires=['scipy', 'natsort', 'pysam', 'python-magic', 'pathos', 'numpy', 'python-levenshtein',
                      'statsmodels']
)