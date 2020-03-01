from setuptools import setup

setup(name='pybustools',
      version=0.1,
      description='Python functions to read the kallisto-bustools format',
      url='http://github.com/redst4r/pybustools/',
      author='redst4r',
      maintainer='redst4r',
      maintainer_email='redst4r@web.de',
      license='GNU GPL 3',
      keywords='RNAseq, kallisto, bustools',
      packages=['pybustools'],
      install_requires=[
          'gmpy2'
          ],
      zip_safe=False)
