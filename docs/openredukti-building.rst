====================
Building OpenRedukti
====================

Build Instructions for Windows
==============================

These are instructions for building OpenRedukti on Windows 10 64-bit.

Pre-Requisites
--------------

OpenRedukti makes use of following external libraries:

* `Protocol Buffers <https://developers.google.com/protocol-buffers/>`_ is used to implement data types
* `OpenBLAS <http://www.openblas.net/>`_ is used for Linear Algebra
* `LAPACK <http://www.netlib.org/lapack/>`_ is used for Linear Algebra
* `CMake <https://cmake.org/>`_ is used to generate build scripts 

Setup OpenBLAS and LAPACK
-------------------------

These are available as pre-built packages from `Ravi Dist <https://github.com/dibyendumajumdar/ravi-dist>`_. 
We assume here that the installed libraries are under ``c:\ravi``.

If you have your OpenBLAS and LAPACK files installed differently, please review and amend the ``FindOpenBLAS.cmake`` file in the ``cmake`` folder.

Build Protocol Buffers
----------------------
I extracted protocol-buffers 3.2 release under ``c:\d\github\protobuf-3.2.0``. 

Since on Windows the release and debug builds depend upon different DLLs, you need to create both versions of protobuf.

I used following steps to create a release build of protobuf::

	cd \d\github\protobuf-3.2.0\cmake
	mkdir buildrelease
	cd buildrelease
	cmake -G "Visual Studio 15 2017 Win64" -DCMAKE_INSTALL_PREFIX=c:\d\protobuf32 -Dprotobuf_MSVC_STATIC_RUNTIME=OFF ..

I then performed Release build using VS2017 followed by INSTALL.

For the debug build I used following steps::

	cd \d\github\protobuf-3.2.0\cmake
	mkdir builddebug
	cd builddebug
	cmake -G "Visual Studio 15 2017 Win64" -DCMAKE_INSTALL_PREFIX=c:\d\protobuf32d -Dprotobuf_MSVC_STATIC_RUNTIME=OFF ..

I then performed a Debug build using VS2017 followed by INSTALL.

After these steps I ended up with following:

* ``c:\d\protobuf32`` - release build
* ``c:\d\protobuf32d`` - debug build

If your installation of Protocol Buffers is different, please review and amend the ``FindPROTOBUF.cmake`` script in the ``cmake`` folder.

Generate Protobuf sources
-------------------------
The next step is to generate the source files for the ``.proto`` definitions. 
From the folder containing OpenRedukti, do following::

	mkdir generated
	mkdir generated\cpp
	cd proto
	generate_protos.bat
	cd ..

This should result in C++ header and source files being generated in ``generated\cpp`` folder.

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
	cmake -DCMAKE_INSTALL_PREFIX=c:\OpenRedukti -G "Visual Studio 15 2017 Win64" -DCMAKE_BUILD_TYPE=Release ..

Remember to select Release configuration in VS2017. You can run the INSTALL target to copy the final binaries to the installation location specified with ``-DCMAKE_INSTALL_PREFIX`` option.

