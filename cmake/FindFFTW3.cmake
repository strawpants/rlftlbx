#CMake module to find the fftw3 libarary

#find fftw3
find_library(FFTW3_LIBRARY NAMES fftw3)

if(${FFTW3_LIBRARY} MATCHES "fftw3")
    set(FFTW3_FOUND TRUE)
endif(${FFTW3_LIBRARY} MATCHES "fftw3")

if(FFTW3_FIND_REQUIRED AND NOT FFTW3_FOUND)
	message( SEND_ERROR "FFTW3 library not found")
endif(FFTW3_FIND_REQUIRED AND NOT FFTW3_FOUND)
