Ncpol2sdpa-Cpp
==
Ncpol2sdpa-Cpp is a C++ library to convert a noncommutative polynomial optimization problem to a sparse semidefinite programming (SDP) problem that can be processed by the [SDPA](http://sdpa.sourceforge.net/) family of solvers. The optimization problem can be unconstrained or constrained by equalities and inequalities.

The objective is to be able to solve very large scale optimization problems. For example, a convergent series of lower bounds can be obtained for ground state problems with arbitrary Hamiltonians.

The implementation has an intuitive syntax for entering Hamiltonians and it scales for a larger number of noncommutative variables using a sparse representation of the SDP problem. 

Dependencies
==
The code requires [SymbolicC++](http://issc.uj.ac.za/symbolic/symbolic.html) to compile and it relies on the C++11 standard. GCC 4.8.1 is known to compile the code. If Ncpol2sdpa-Cpp is compiled with OpenMP support, SymbolicC++ needs a [patch](http://peterwittek.com/files/openmp_patch.txt) to ensure thread-safety.

Usage
==
A simple usage example is included in examplencpol.cpp. A more sophisticated application is given in hamiltonian.cpp, which implements the Hamiltonian of a fermionic system in a 2D grid.

The implementation installs as a library. Subsequent use must specify the include directory of the header files and the library for compilation. 

Compilation & Installation
==
From GIT repository first run

    $ ./autogen.sh

Then follow the standard procedure:

    $ ./configure [options]
    $ make
    $ make install

Options for configure

    --enable-openmp Enable OpenMP support (experimental)

OpenMP support is still experimental and deadlocks occur in larger problems.

    --with-symbolicc++-incdir=DIR   SymbolicC++ include directory [default /usr/include]
    --with-symbolicc++-libdir=DIR   SymbolicC++ library directory [default /usr/lib]

Specify these directories if the compiler cannot see them in the include and library paths.

Known Issues
==
Hermicity of noncommuting variables is not handled correctly.

The fast submonomial substitution does not handle some rare cases.

More Information
==
For more information refer to the following manuscript:

[http://arxiv.org/abs/1308.6029](http://arxiv.org/abs/1308.6029)
