try:
    from setuptools import setup
except:
    from distutils.core import setup

setup(name = 'bpch',
      version = '1.0rc',
      author = 'Barron Henderson',
      author_email = 'barronh@ufl.edu',
      maintainer = 'Barron Henderson',
      maintainer_email = 'barronh@ufl.edu',
      py_modules = ['bpch'],
      install_requires = ['numpy (>=1.5)', 'matplotlib (>=1.0)'],
      requires = ['numpy (>=1.5)', 'matplotlib (>=1.0)'],
      url = 'https://github.com/barronh/bpch',
      download_url = 'https://github.com/barronh/bpch/archive/1.0rc.zip'
      )
