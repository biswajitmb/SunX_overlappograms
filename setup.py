import os
from setuptools import setup, find_packages

setup(name='SunX_overlappograms',
      version='0.1',
      description='SunX_overlapograms: A python package to simulate the overlappograms from emission measure cube',
      author='Biswajit Mondal',
      author_email='biswajit70mondal94@gmail.com',
      url='https://github.com/biswajitmb/SunX_overlappograms.git',
      packages=find_packages(),
      install_requires=[
        'numpy',
        'astropy',
        'scipy',
        'matplotlib',
        'sunpy',
        'aiapy',
      ],
     )
