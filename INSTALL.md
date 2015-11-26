# HOPR INSTALLATION PROCEDURE


## Prerequisites

HOPR supports Linux-based systems only requires a x86_64
compliant platform and has been tested on the following platforms:

- Ubuntu 12.04 or newer
- Linux Mint 17 or newer
- SUSE Linux Enterprise Server 11 SP3 and newer


### Compilers

HOPR requires a C and a Fortran 2003 compliant compiler,
compilers tested with HOPR include

- GNU Compiler Collection 4.6 or newer
- Intel C/Fortran Compiler 12 or newer (recommended)
- CRAY Compiler Environment 8.1 or newer


### Libraries

The following libraries are required, if not mentioned
otherwise, including their development headers:

- libc6
- zlib
- BLAS/LAPACK (or compatible, e.g. ATLAS, MKL)
- Python 2.7 or newer (optional)

If not present on your system, HOPR will automatically
download and compile these libraries

- HDF5
- CGNS


## Compiling HOPR

To compile HOPR, specify the compiler in the file `Makefile.defs`,

- Intel: `COMPILER=INTEL`
- GNU:   `COMPILER=GNU`

then compile HOPR from the command-line with `make`. This will
also download and compile all libraries in the share folder.

To compile HOPR with MPI support (e.g. to link other libraries
using MPI) or with debug flags add the Flag `_MPI` or `_DEBUG`
to the `COMPILER` string, e.g. `COMPILER=GNU_MPI_DEBUG`.
Other options are descibed in the `Makefile.defs` header.

In case you want to use a precompiled HDF5 version on your system,
set the option `BUILD_HDF5=<empty>`, then the path to HDF5 can be
specified by the environment variable `$HDF5_DIR`.


## Testing HOPR

After compiling, you can test HOPR by going to the
`tutorials` directory and running the script `executeall.sh`,
which will run HOPR for all tutorial cases. For further
information check the [HOPR website](http://www.hopr-project.org/).

