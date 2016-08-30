#CMake module to find the fortran nectdf interface

#find netcdf fortran interface
find_library(NETCDFF_LIBRARY NAMES netcdff)

if(${NETCDFF_LIBRARY} MATCHES "netcdff")
    set(NETCDFF_FOUND TRUE)
    get_filename_component(NETCDFF_ROOTPATH "${NETCDFF_LIBRARY}" DIRECTORY)
    set(NETCDFF_MODPATH ${NETCDFF_ROOTPATH}/../include CACHE PATH "netcdf module path")
endif(${NETCDFF_LIBRARY} MATCHES "netcdff")

if(NETCDFF_FIND_REQUIRED AND NOT NETCDFF_FOUND)
	message( SEND_ERROR "Netcdf Fortran interface not found")
endif(NETCDFF_FIND_REQUIRED AND NOT NETCDFF_FOUND)
