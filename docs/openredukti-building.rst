====================
Building OpenRedukti
====================

Dependencies
------------

OpenRedukti makes use of following external libraries:

* `Protocol Buffers <https://developers.google.com/protocol-buffers/>`_ is used to implement data types
* `OpenBLAS <http://www.openblas.net/>`_ is used for Linear Algebra
* `LAPACK <http://www.netlib.org/lapack/>`_ is used for Linear Algebra
* `CMake <https://cmake.org/>`_ is used to generate build scripts 

Build Instructions for Ubuntu Linux 18.04 LTS
=============================================

Pre-Requisites
--------------

Install following::

    sudo apt install git
    sudo apt install cmake
    sudo apt install libreadline-dev
    sudo apt install libopenblas-dev
    sudo apt install libprotobuf-dev
    sudo apt install protobuf-compiler

Build OpenRedukti
-----------------

Clone the OpenRedukti github repository and do following:: 

    mkdir buildrelease
    cd buildrelease
    cmake -DCMAKE_INSTALL_PREFIX=~/Software/OpenRedukti -DCMAKE_BUILD_TYPE=Release ..
    make install

Build Instructions for Windows
==============================
These are instructions for building OpenRedukti on Windows 10 64-bit.


Setup OpenBLAS and LAPACK
-------------------------
These are available as pre-built packages from `Ravi Distribution Dependencies <https://github.com/dibyendumajumdar/ravi-external-libs>`_. 
We assume here that the installed libraries are under ``c:\Software\OpenRedukti``. 

If you have your OpenBLAS and LAPACK files installed differently, please review and amend the ``FindOpenBLAS.cmake`` file in the ``cmake`` folder.

Obtain Protocol Buffers via vcpkg
---------------------------------
Install `vcpkg <https://github.com/Microsoft/vcpkg>`_.
We assume below that ``vcpkg`` is installed at ``c:\work\vcpkg``.

Get protobuf as follows::

    vcpkg install protobuf:x64-windows

On my machine after installation I get this::

    C:\work\vcpkg>vcpkg list
    protobuf:x64-windows                               3.6.1-2          Protocol Buffers - Google's data interchange format

Ensure protoc is on the path as follows::

    set PATH=C:\work\vcpkg\installed\x64-windows\tools\protobuf;%PATH%


Build OpenRedukti
-----------------
Once all of above steps are done, you can build OpenRedukti as follows::

	mkdir build
	cd build
	cmake -G "Visual Studio 15 2017 Win64" -DCMAKE_BUILD_TYPE=Debug ..

Above creates projects suited for debug build. You can go into VS2017 and do the build from there.

For a release build, do following::

	mkdir buildrelease
	cd buildrelease
	cmake -DCMAKE_INSTALL_PREFIX=c:\Software\OpenRedukti -G "Visual Studio 15 2017 Win64" -DCMAKE_BUILD_TYPE=Release ..

Remember to select Release configuration in VS2017. You can run the INSTALL target to copy the final binaries to the installation location specified with ``-DCMAKE_INSTALL_PREFIX`` option.

