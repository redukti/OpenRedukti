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

Optionally if you want to enable a gRPC based server application then additional dependency on:

* `gRPC <https://grpc.io/>`_ framework

Note that the gRPC server is required if you want to use the pricing and curve building functions from the `Python interface <https://github.com/redukti/PyRedukti>`_.

Build Instructions for RHEL 7.6
===============================
I am using gcc 8.3 on Redhat obtained via devtoolset-8.
I installed CMake manually.

I built protobuf and grpc manually (see below).

I enabled `EPEL <https://fedoraproject.org/wiki/EPEL>`_ repository. This was necessary to obtain openblas which I then installed as follows::

  sudo yum install openblas-devel

Building protobuf
-----------------

After downloading and unpacking the sources I executed::

  autogen.sh
  configure --prefix=$HOME/Software/protobuf
  make
  make install
	
Building GRPC
-------------

I tried building GRPC from the release packages but this failed due to unsatisfied dependency on libcares. Somehow the default method is not compatible with RHEL.

So then I cloned the `gRPC github repo <https://github.com/grpc/grpc>`_ and built from there.

Note that I had to modify following changes:

In the Supplied ``Makefile``:

* ``prefix`` - base path for installation
* ``CPPFLAGS`` - I removed the ``-Werror`` option as this caused a failure with gcc 8.3

In ``third_party/boringssl/CMakeLists.txt``:

* ``C_CXX_FLAGS`` - I removed the ``-Werror`` option as this caused a failure with gcc 8.3

I executed following steps::

  git submodule update --init
  # Make changes above
  make
  make install
	
The installation churned out couple of permission denied messages but succeeded.

Building OpenRedukti
--------------------

OpenRedukti was built as follows::

  mkdir build
  cd build
  cmake -DGRPC_SERVER=ON -DProtobuf_ROOT=~/Software/protobuf -DCMAKE_INSTALL_PREFIX=~/Software/redukti -DGRPC_ROOT=~/Software/grpc ..

You can omit the ``GRPC_SERVER`` option if that is not needed.

Build Instructions for Ubuntu Linux 18.04 LTS
=============================================

Following instructions do not build GRPC server.

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
NOTE: I had to build protobuf locally because I faced some issues with below. (FIXME)

Install `vcpkg <https://github.com/Microsoft/vcpkg>`_.
We assume below that ``vcpkg`` is installed at ``c:\work\vcpkg``.

Get protobuf as follows::

    vcpkg install protobuf:x64-windows

On my machine after installation I get this::

    C:\work\vcpkg>vcpkg list
    protobuf:x64-windows                               3.6.1-2          Protocol Buffers - Google's data interchange format

Ensure protoc is on the path as follows::

    set PATH=C:\work\vcpkg\installed\x64-windows\tools\protobuf;%PATH%

Build gPRC
----------
This is an optional step. 

On Windows, you can build and install gRPC using `vcpkg`. This is what I did.
Or else follow instructions at `gRPC C++ Building from source <https://github.com/grpc/grpc/blob/master/BUILDING.md>`_.  

Build OpenRedukti
-----------------
Once all of above steps are done, you can build OpenRedukti as follows::

	mkdir build
	cd build
	set PATH=c:\Software\protobuf371d\bin;%PATH%
	cmake  -DCMAKE_INSTALL_PREFIX=c:\Software\OpenRedukti -G "Visual Studio 15 2017 Win64" -DCMAKE_BUILD_TYPE=Debug -DPROTOBUF_SRC_ROOT_FOLDER=c:\Software\protobuf371d -DgRPC_DIR=c:\work\vcpkg\installed\x64-windows-static-dyncrt\share\grpc -Dc-ares_DIR=c:\work\vcpkg\installed\x64-windows-static-dyncrt\share\c-ares ..

Above creates projects suited for debug build. You can go into VS2017 and do the build from there.

For a release build, do following::

	mkdir buildrelease
	cd buildrelease
	cmake -DCMAKE_INSTALL_PREFIX=c:\Software\OpenRedukti -G "Visual Studio 15 2017 Win64" -DCMAKE_BUILD_TYPE=Release ..

Remember to select Release configuration in VS2017. You can run the INSTALL target to copy the final binaries to the installation location specified with ``-DCMAKE_INSTALL_PREFIX`` option.

