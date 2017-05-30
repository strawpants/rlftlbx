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

### Clone rlftbx
```
git@gitlab.igg.uni-bonn.de:roelof/rlftlbx.git

```

### Create a build directory and run cmake 
```
mkdir rlftlbx_build;cd rlftblx_build

cmake -DCMAKE_PREFIX_PATH="/

```

