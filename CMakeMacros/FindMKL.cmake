
if (MKL_LIBRARIES)
  set(MKL_FIND_QUIETLY TRUE)
endif (MKL_LIBRARIES)

if(CMAKE_MINOR_VERSION GREATER 4)

if(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")

find_library(MKL_LIBRARIES
  mkl_core
  PATHS
  $ENV{MKLLIB}
  /opt/intel/*/mkl/lib/intel64
  ${LIB_INSTALL_DIR}
)

if(MKL_LIBRARIES)
  set(MKL_LIBRARIES ${MKL_LIBRARIES} mkl_intel_lp64 mkl_core mkl_intel_thread iomp5 dl pthread m ) # mkl_gnu_thread dl  or  mkl_intel_thread iomp5 dl  or  mkl_sequential 
endif()

else(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")

find_library(MKL_LIBRARIES
  mkl_core
  PATHS
  $ENV{MKLLIB}
  /opt/intel/mkl/*/lib/32
  /opt/intel/Compiler/*/*/mkl/lib/32
  ${LIB_INSTALL_DIR}
)

find_library(MKL_GUIDE
  guide
  PATHS
  $ENV{MKLLIB}
  /opt/intel/mkl/*/lib/32
  /opt/intel/Compiler/*/*/mkl/lib/32
  /opt/intel/Compiler/*/*/lib/intel32
  ${LIB_INSTALL_DIR}
)

if(MKL_LIBRARIES)
  set(MKL_LIBRARIES ${MKL_LIBRARIES} mkl_intel mkl_sequential ${MKL_GUIDE} pthread)
endif()

endif(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")

endif(CMAKE_MINOR_VERSION GREATER 4)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MKL DEFAULT_MSG MKL_LIBRARIES)

mark_as_advanced(MKL_LIBRARIES)
