from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='modelbase',
      version='0.2.3',
      description='A package to build metabolic models',
      long_description=long_description,
      url='https://gitlab.com/ebenhoeh/modelbase',
      author='Oliver Ebenhoeh',
      author_email='oliver.ebenhoeh@hhu.de',
      license='GPL4',
      packages=['modelbase'],
      install_requires=[
          'numpy',
          'scipy',
          'matplotlib'
          #'itertools',
          #'re',
          #'pickle'
      ],
      zip_safe=False)
