#CMake module to find an optimzed BLAS Library (the default Module FindBLAS.cmake has some problems and can't find Openblas
#this also sets lapack libraries if they are part of the library
#not to find openblas or atlas in non-standard directories pass -D CMAKE_PREFIX_PATH=/weird/directory to cmake

list(APPEND CMAKE_FIND_LIBRARY_SUFFIXES .so.0)
find_library(OPTBLAS_LIBRARY NAMES openblas atlas)
set(BLAS_libs "")

if(${OPTBLAS_LIBRARY} MATCHES "atlas")
    set(OPTBLAS_VENDOR "ATLAS")
    
    #also try to find lapack (usually included)
    find_library(OPTLAPACK_LIBRARY NAMES lapack)
    if(OPTLAPACK_LIBRARY)
        set(OPTLAPACK_FOUND TRUE)
    list(APPEND BLAS_libs "-llapack")
    endif(OPTLAPACK_LIBRARY)
    list(APPEND BLAS_libs "-llapack" "-lptf77blas" "-lptcblas" "-latlas" "-lpthread")
elseif(${OPTBLAS_LIBRARY} MATCHES "openblas")
    set(OPTBLAS_VENDOR "OpenBLAS")
    list(APPEND BLAS_libs ${OPTBLAS_LIBRARY})
    set(OPTLAPACK_FOUND TRUE)
endif(${OPTBLAS_LIBRARY} MATCHES "atlas")

set( OPTBLAS_LIBRARIES ${BLAS_libs} CACHE STRING "BLAS and Lapack librarries to include")
get_filename_component(OPTBLAS_LIBPATH ${OPTBLAS_LIBRARY} DIRECTORY CACHE)

if(OPTBLAS_VENDOR)
    message("Optimized BLAS found: " ${OPTBLAS_VENDOR} " " ${OPTBLAS_LIBRARIES} )
    set(OPTBLAS_FOUND TRUE)
endif(OPTBLAS_VENDOR)
#little clean up of CMakeCache.txt
unset(OPTBLAS_LIBRARY CACHE)

