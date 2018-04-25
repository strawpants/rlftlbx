## Fortran command line tools to do  geodetic stuff.

This software suite is my heritage from my PhD time. It is a set of classical  
command line tools based upon a small library, which do relatively small tasks.  
These are stitched together with for example bash scripts when more complex  
problems are to be solved.  

Although I still do a significant amount of work with these tools, I intend to  
gradually replace it with Frommmle, a c++ library suite which can be scaled more  
easily to parallel computation environments, and which is (hoperfully) better  
maintainable.  


## Installation

prerequisites:

* install netcdf-fortran (and friends: hdf5,netcdf)
* install [shtools](https://github.com/SHTOOLS/SHTOOLS)
* install lapack/blas library

### A note on building/linking netcdf-fortran
It's highly advisable to create a central place where you install both the netcdf-fortran *and* the netcdf-c libraries to. So this means either installing both netcdf interfaces to the system default `/usr/lib` `/usr/include` or in a custom folder (e.g. `/opt/roelof/lib` and `/opt/roelof/include`). Problems will start to occur when:
* the C library (libnetcdf.so) and the fortran library live in separate directories
* include files (netcdf.h netcdf.mod) live in separate directories
* The system netcdf-c library is in conflict with the needed c library required by the fortran interface


### Clone rlftbx
```
git@gitlab.igg.uni-bonn.de:roelof/rlftlbx.git

```

### Create a build directory and run cmake 
```
mkdir rlftlbx_build;cd rlftblx_build

cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="ADDITIONALLIBPATH1;PATH2" -DCMAKE_INSTALL_PREFIX="INSTALLPATH" PATHTORLFTLBXROOT
```
Alternatively if this does not work try to be more specific:
```
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_LIBRARY_PATH="/path/to/BLAS/lib;/path/to/netcdf/fortran/lib64;/path/to/SHTOOLS/lib" -DCMAKE_INCLUDE_PATH="/path/to/netcdf/include;/path/to/BLAS/include" -DCMAKE_INSTALL_PREFIX="INSTALLPATH" PATHTORLFTLBXROOT
```



Once cmake completes successfully, run:
```
make 

make install

```

 



