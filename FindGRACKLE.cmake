# - Find the GRACKLE library
#
# Usage:
#   find_package(GRACKLE [REQUIRED] [QUIET] )
#
# It sets the following variables:
#   GRACKLE_FOUND               ... true if grackle is found on the system
#   GRACKLE_LIBRARY             ... full path to grackle library
#   GRACKLE_INCLUDE             ... grackle include directory
#
# The following variables will be checked by the function
#   GRACKLE_ROOT               ... path to the GRACKLE install directory,
#                              ... usually the MACH_INSTALL_PREFIX when installing it
#

#If environment variable GRACKLE_DIR is specified, it has same effect as GRACKLE_ROOT
if( NOT GRACKLE_ROOT AND ENV{GRACKLE_DIR} )
  set( GRACKLE_ROOT $ENV{GRACKLEDIR} )
endif()

# Check if we can use PkgConfig
find_package(PkgConfig)

#Determine from PKG
if( PKG_CONFIG_FOUND AND NOT GRACKLE_ROOT )
  pkg_check_modules( PKG_GRACKLE QUIET "grackle" )
endif()


if( GRACKLE_ROOT )

  find_library(
    GRACKLE_LIB
    NAMES "grackle"
    PATHS ${GRACKLE_ROOT}
    PATH_SUFFIXES "lib"
    NO_DEFAULT_PATH
  )

  find_path(
    GRACKLE_INCLUDE
    NAMES "grackle.h"
    PATHS ${GRACKLE_ROOT}
    PATH_SUFFIXES "include"
    NO_DEFAULT_PATH
  )

endif( GRACKLE_ROOT )

set(GRACKLE_LIBRARY ${GRACKLE_LIB})


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GRACKLE DEFAULT_MSG
                                  GRACKLE_INCLUDE GRACKLE_LIBRARY)

mark_as_advanced(GRACKLE_INCLUDE GRACKLE_LIBRARY GRACKLE_LIB)
