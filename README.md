## Fortran command line tools to do  geodetic stuff.

This software suite orginates from my PhD time. It is a set of classical unix-like 
command line tools based upon a small library, which do relatively small tasks.  
These are stitched together with for example bash scripts when more complex problems are to be solved.

Although I still do some work with these tools, I intend to  
gradually replace it with [Frommle](https://github.com/strawpants/frommle), a c++ library suite with a python facing interface. One of the reasons is that, even though functionaility has been put into libraries, some programs have still a rather monolithic structure, and would be tedious to refactor. Futhermore, I found myself being limited by the fortran language toolset, and lost interest in actively developing tools with it.

TLDR; 
This repository will receive small bugfixes when the need arises, but I don't intent to actively develop it anymore.


## Installation

prerequisites:

* install netcdf-fortran (and friends: hdf5,netcdf)
* install [shtools](https://github.com/SHTOOLS/SHTOOLS)
* install a capable lapack/blas library such as OpenBLAS

### A note on building/linking netcdf-fortran
It's highly advisable to create a central place where you install both the netcdf-fortran *and* the netcdf-c libraries to. So this means either installing both netcdf interfaces to the system default `/usr/lib` `/usr/include` or in a custom folder (e.g. `/opt/roelof/lib` and `/opt/roelof/include`). Problems will start to occur when:
* the C library (libnetcdf.so) and the fortran library live in separate directories
* include files (netcdf.h netcdf.mod) live in separate directories
* The system netcdf-c library is in conflict with the needed c library required by the fortran interface


### Clone rlftbx
```
 git clone git@github.com:strawpants/rlftlbx.git

```

### Create a build directory (preferabbly outside of the source tree) and run cmake 
```
mkdir rlftlbx_build;cd rlftblx_build

cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="ADDITIONALLIBPATH1;PATH2" -DCMAKE_INSTALL_PREFIX="INSTALLPATH" PATHTORLFTLBXROOT
```
Alternatively if this doesi, and libaries are not found one can try to give cmake more hints on where to find the libs:
```
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_LIBRARY_PATH="/path/to/BLAS/lib;/path/to/netcdf/fortran/lib64;/path/to/SHTOOLS/lib" -DCMAKE_INCLUDE_PATH="/path/to/netcdf/include;/path/to/BLAS/include" -DCMAKE_INSTALL_PREFIX="INSTALLPATH" PATHTORLFTLBXROOT
```



Once cmake completes successfully, run:
```
make 
make install

```

 



