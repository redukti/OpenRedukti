============
Introduction
============

OpenRedukti is a C++ library for working with Interest Rate Derivative products such as Interest Rate Swaps, and
FRAs. It allows you to build Interest Rate curves with different interpolation methods, and then use these curves
to compute present value and sensitivities of Interest Rate Derivatives.

OpenRedukti is Free Software, licensed under the GNU General Public License, v3. 

Main Features
=============
* Small library with minimal external dependencies (only external dependencies are BLAS, LAPACK and Google Protocol Buffers) 
* Ability to express an interest rate product as a set of cashflows
* Bootstrap continuously compounded zero coupon interest rate curves using Linear, CubicSpline, and MonotoneConvex interpolators
* Interpolate curves in the discount factor space using LogLinear and LogCubicSpline interpolators
* Compute present value of cashflows
* Compute first and second order derivatives using `Automatic/algorithmic Differentiation <http://www.autodiff.org/>`_.
* Script using `Ravi <https://github.com/dibyendumajumdar/ravi>`_ - a derivative of `Lua <http://www.lua.org>`_ programming language
* **New!** Script using Python - see `PyRedukti <https://github.com/redukti/PyRedukti>`_ for details.
* **New!** `gRPC <https://grpc.io/>`_ server implementation enables access from multiple programming languages.


Background
==========
OpenRedukti was part of the MyCCP product that was being developed by REDUKTI LIMITED. Sadly development of MyCCP ceased in 2017 due to lack of funding.

The main differences between the Open Source release and the proprietary version used in MyCCP are:

* Only the core C++ pricing library has been released
* The functionality for generating cashflows from FpML trades has not been released as this is fine tuned for the needs of MyCCP
* The Limit Checker and VaR Calculator have not been released
* The MyCCP front-end and middle tier components, written in C#, have not been released as these are very specific to requirements of a CCP.

Documentation
=============

See `redukti.github.io <https://redukti.github.io/>`_. 

Ackowledgements
===============
OpenRedukti gratefully acknowledges ideas and code it is using from other projects.

* My good friend `Christer Rydberg <https://www.linkedin.com/in/christer-rydberg-phd-98012a7/>`_ for showing me `how to implement Automatic Differentiation <https://github.com/redukti/OpenRedukti/blob/master/docs/Sensitivities.pdf>`_ and for help and advice over the years. 
* `Jeffrey Fike's work <http://adl.stanford.edu/hyperdual/>`_ on automatic differentiation using hyperdual vectors.
* `QuantLib <http://quantlib.org/index.shtml>`_ which has provided the basis for some of the key components such as interpolators.
* The monotone convex interpolation method is based on papers and VB code from `Financial Modeling Agency <http://finmod.co.za/#our-research>`_. 
* The `C/C++ Minpack <http://devernay.free.fr/hacks/cminpack/>`_ library provides the Levenberg-Marquardt solver used in curve building.
