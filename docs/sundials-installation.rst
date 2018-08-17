Building Sundials
=================

Prerequisites:
* Ubuntu: gcc, build-essential, cmake, cmake-curses-gui

Download sundials v2.6.0
.. code-block:: unix
   cd sundials
   mkdir instdir
   mkdir builddir
   cd builddir
   cmake ../srcdir
   ccmake ../srcdir

In ccmake set

.. code-block:: unix
   CMAKE_C_FLAGS=-fPIC
   LAPACK_ENABLE=OFF
   CMAKE_INSTALL_PREFIX = /home/marvin/sundials/instdir
   EXAMPLES_INSTALL_PATH = =/home/marvin/sundials/instdir/examples

Installing Assimulo
===================
.. code-block:: python3
   conda install cython
   conda install nose
   conda install -c chria assimulo
