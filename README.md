# About
This repository holds a compilable version of the [Mitsuba Renderer](https://github.com/mitsuba-renderer/mitsuba). The original code and build scripts are somewhat outdated, and do not work quite well on today's machine.

# Note
Currently, only Linux system is supported.

# How to Compile
```shell
# Install utilities
sudo apt install git cmake ninja-build clang

# Install dependencies
sudo apt install libxerces-c-dev libboost-all-dev libeigen3-dev libglewmx-dev libxxf86vm-dev libpng-dev libjpeg-dev libopenexr-dev libcollada-dom-dev libfftw3-dev qtbase5-dev qtchooser qt5-qmake qtbase5-dev-tools libqt5xmlpatterns5-dev libpython-all-dev python3 libpython3-all-dev

# Clone the repository
git clone https://github.com/yysun997/mitsuba.git

# Configure and compile (Configure warning should not happen, but compile warning is fine.)
mkdir ./mitsuba/build/release
cd ./mitsuba/build/release
cmake -G Ninja -DCMAKE_CXX_COMPILER=clang -DCMAKE_C_COMPILER=clang -DCMAKE_BUILD_TYPE=Release ../..
ninja

# Take your python version
ls

# The command-line
# ./mitsuba -h

# The GUI
# ./mtsgui

# The python bindings
# python3
# >>> import sys
# >>> sys.path.append("./python${Python_VERSION}")
# >>> import mitsuba
# >>> from mitsuba.core import *
# >>> vec = Vector(1, 2, 3)
# >>> normalize(vec)

# See official documentation for more details
# https://www.mitsuba-renderer.org/releases/current/documentation.pdf
```
