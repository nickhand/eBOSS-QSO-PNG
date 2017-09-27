from setuptools import setup, find_packages
import os

setup(
    name='eboss_qso',
    version='0.0.1',
    author='Nick Hand',
    packages=find_packages(),
    description=("Python code to analyze primordial non-Gaussianity "
                 "from the eBOSS QSO sample"),
    install_requires=['numpy', 'scipy']
)
