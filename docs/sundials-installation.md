# Building sundials

## Prerequisites


* Ubuntu
  * gcc
  * build-essential
  * cmake
  * cmake-curses-gui

Download sundials v2.6.0

```
  cd sundials
  mkdir instdir
  mkdir builddir
  cd builddir
  cmake ../srcdir
  ccmake ../srcdir
```

In ccmake set

```
  CMAKE_C_FLAGS=-fPIC
  LAPACK_ENABLE=OFF
  CMAKE_INSTALL_PREFIX = /home/marvin/sundials/instdir
  EXAMPLES_INSTALL_PATH = =/home/marvin/sundials/instdir/examples
```

Installing Assimulo
===================
```python3
   conda install cython
   conda install nose
   conda install -c chria assimulo
```
