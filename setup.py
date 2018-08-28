# This file is part of modelbase.
#
# modelbase is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# modelbase is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with modelbase.  If not, see <http://www.gnu.org/licenses/>.

from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='modelbase',
      version='0.2.4',
      description='A package to build metabolic models',
      long_description=long_description,
      url='https://gitlab.com/ebenhoeh/modelbase',
      author='Oliver Ebenhoeh',
      author_email='oliver.ebenhoeh@hhu.de',
      license='GPL3',
      packages=['modelbase'],
      install_requires=[
          'numpy',
          'scipy',
          'matplotlib',
          'pandas',
          'numdifftools'
      ],
      zip_safe=False)
