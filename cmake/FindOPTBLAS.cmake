#CMake module to find an optimzed BLAS Library (the default Module FindBLAS.cmake has some problems and can't find Openblas
#this also sets lapack libraries if they are part of the library
#not to find openblas or atlas in non-standard directories pass -D CMAKE_PREFIX_PATH=/weird/directory to cmake

list(APPEND CMAKE_FIND_LIBRARY_SUFFIXES .so.0)
find_library(OPTBLAS_LIBRARIES NAMES openblas atlas)

if(${OPTBLAS_LIBRARIES} MATCHES "atlas")
    set(OPTBLAS_VENDOR "ATLAS")
    
    #also try to find lapack (usually included)
    find_library(OPTLAPACK_LIBRARIES NAMES lapack)
    if(OPTLAPACK_LIBRARIES)
        set(OPTLAPACK_FOUND TRUE)
    endif(OPTLAPACK_LIBRARIES)
    set(GEOAUX_LIBS ${GEOAUX_LIBS} "-lptf77blas.a -lpthread")
elseif(${OPTBLAS_LIBRARIES} MATCHES "openblas")
    set(OPTBLAS_VENDOR "OpenBLAS")
    set(OPTLAPACK_FOUND TRUE)
endif(${OPTBLAS_LIBRARIES} MATCHES "atlas")
  
if(OPTBLAS_VENDOR)
    message("Optimized BLAS found: " ${OPTBLAS_VENDOR} " " ${OPTBLAS_LIBRARIES})
    set(OPTBLAS_FOUND TRUE)
endif(OPTBLAS_VENDOR)
