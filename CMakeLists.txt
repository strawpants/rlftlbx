cmake_minimum_required(VERSION 2.8.12)
project(RLFtlbx)

#set some options
option(USE_OPENMP "Build with OPENMP enabled" OFF)


enable_language(Fortran)


# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING   "Choose the type of build, options are: None Debug Release." FORCE)
endif (NOT CMAKE_BUILD_TYPE)



#basic compiler flags
set(CMAKE_Fortran_FLAGS_RELEASE "-march=native -O3")
set(CMAKE_Fortran_FLAGS_DEBUG "-g") 

##add some flags and libraries when using openmp
if( ${USE_OPENMP} )
    find_package(OpenMP REQUIRED)
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${OpenMP_Fortran_FLAGS}" )
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${OpenMP_Fortran_FLAGS}" )
endif( ${USE_OPENMP} )


#add make info
set(CMAKE_VERBOSE_MAKEFILE on)


#set the name of the library to be build
set(RLFLIB "RLFtlbx")
set(RLFLIBMAJORVERSION 0)
set(RLFLIBMINORVERSION 1)

#where to put build executables
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

#where to put archives and the libraries
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

#find prerequisites (BLAS & LAPACK)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
#try finding an optimizedBLAS library
find_package(OPTBLAS REQUIRED)

#find the netcdf library
find_package(NetCDFF REQUIRED )

#find the SHTOOLS library
find_package(SHTOOLS REQUIRED )

#find the fftw3 libarary
find_package(FFTW3 REQUIRED )

include_directories("${SHTOOLS_MODPATH}" "${NETCDFF_MODPATH}")



#add module sources
LIST(APPEND F90MODULES FORTtlbx/AutoRegressive.f90 FORTtlbx/BINfiletools.f90 FORTtlbx/EOP_tools.f90 FORTtlbx/FORTtlbx.f90 FORTtlbx/GroopsFile.f90 FORTtlbx/BIN_operations.f90 )
LIST(APPEND F90MODULES GPStlbx/GPStlbx.f90 GPStlbx/SINEXtools.f90 )
LIST(APPEND F90MODULES GRDtlbx/grdtlbx.f90 GRDtlbx/PSMSL_TOOL.f90 )
LIST(APPEND F90MODULES SHtlbx/SHtlbx.f90 SHtlbx/SH_trpnt_mem.f90 SHtlbx/SH_synthesis.f90 SHtlbx/GIA_func.f90 )
add_library(mod OBJECT ${F90MODULES} ${COREHEADERS})

#add additional library source files
LIST(APPEND F90OBJS FORTtlbx/BIN_MATH.f90 FORTtlbx/BIN_MAT_TIME.f90 FORTtlbx/BLAS_custom.f90 FORTtlbx/GRACENORMtools.f90 FORTtlbx/NORM_fix_red.f90 FORTtlbx/regexf.f90 FORTtlbx/Sort.f90 FORTtlbx/statistics.f90 FORTtlbx/various.f90)

LIST(APPEND F90OBJS GPStlbx/GPS_XYZ_geo_trans.f90 GPStlbx/get_sinex.f90 GPStlbx/GPS_rotmat.f90 GPStlbx/get_sinexblock.f90 GPStlbx/GPS_time.f90 GPStlbx/write_sinex.f90 GPStlbx/GPS_discont_cumu.f90 GPStlbx/conv2itrf2000.f90 GPStlbx/GPS_obseq.f90)

LIST(APPEND F90OBJS SHtlbx/rotfeedback.f90 SHtlbx/index.f90 SHtlbx/SH_gaus.f90 SHtlbx/SH_read.f90 SHtlbx/SH_geocmod.f90 SHtlbx/SH_write.f90 SHtlbx/SH_loadlove.f90 SHtlbx/SH_tripint.f90  SHtlbx/SH_overlapint.f90 SHtlbx/c20_replace.f90 SHtlbx/SH_loadGIA_uplift_ratio.f90 SHtlbx/SH_Love_d1_trans.f90 )

LIST (APPEND F90OBJS GRDtlbx/NC_write.f90 GRDtlbx/OBP_obseq.f90)


#Cstuff 
LIST(APPEND COBJS FORTtlbx/get_permvecC.c  FORTtlbx/REGEX.c  FORTtlbx/Sort_various.c  FORTtlbx/STREAMFILE.c )


add_library(${RLFLIB} STATIC $<TARGET_OBJECTS:mod> ${F90OBJS} ${COBJS} )
install(TARGETS ${RLFLIB} DESTINATION lib)


#add programs
set(FORTPROG ascii_2_BIN BIN_BDmerge BIN_BDsplit BIN_buildnorm BIN_MATH BIN_MAT_TIME BIN_swiss fit1D GRACE2BIN Groops2BIN NORM_tool)
foreach(prog  ${FORTPROG})
add_executable(${prog} FORTtlbx/${prog}.f90)
target_link_libraries(${prog} ${RLFLIB} ${OPTBLAS_LIBRARIES})
install(TARGETS ${prog} DESTINATION bin)
endforeach(prog in ${FORTPROG})


set(GPSPROG SNX2BIN GPS_calendar GPS_designmat)
foreach(prog  ${GPSPROG})
add_executable(${prog} GPStlbx/${prog}.f90)
target_link_libraries(${prog} ${RLFLIB} ${OPTBLAS_LIBRARIES} ${SHTOOLS_LIBRARY})
install(TARGETS ${prog} DESTINATION bin)
endforeach(prog in ${GPSPROG})

set(SHPROG SH_fit JKcov_2_BIN SH_decormatSW SH_prod2sum_mat_openmp SH_transmat SHiso_2_colat SH_axisload SH_designmat SH_sealevmat nc_2_SH SH_2_nc SH_dump SH_basinav SH_filter SH_prod2sum_mat SH_diagreg SH_degvar SH_2_pnt SHSELEN_xtract SH_gradMAT xyz_2_SH BIN_2_JKmat)
foreach(prog  ${SHPROG})
add_executable(${prog} SHtlbx/${prog}.f90)
target_link_libraries(${prog} ${RLFLIB} ${OPTBLAS_LIBRARIES} ${SHTOOLS_LIBRARY} ${NETCDFF_LIBRARY} ${FFTW3_LIBRARY})
install(TARGETS ${prog} DESTINATION bin)
endforeach(prog in ${SHPROG})


set(GRDPROG OBP_normeq PSMSL_READ )
foreach(prog  ${GRDPROG})
add_executable(${prog} GRDtlbx/${prog}.f90)
target_link_libraries(${prog} ${RLFLIB} ${OPTBLAS_LIBRARIES} ${NETCDFF_LIBRARY} ${SHTOOLS_LIBRARY})
install(TARGETS ${prog} DESTINATION bin)
endforeach(prog in ${GRDPROG})


