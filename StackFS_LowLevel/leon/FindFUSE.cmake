# This module can find FUSE Library
#
# The following variables will be defined for your use:
# - FUSE_FOUND : was FUSE found?
# - FUSE_INCLUDE_DIRS : FUSE include directory
# - FUSE_LIBRARIES : FUSE library

find_path(
    FUSE_INCLUDE_DIRS
    NAMES fuse_common.h fuse_lowlevel.h fuse.h
    PATHS /var/services/homes/gbalakrishna/fuse3/fuse-library/include /usr/lib/fuse-library/include)

set(fuse_names fuse3)
if(APPLE)
  list(APPEND fuse_names libosxfuse.dylib)
endif()

find_library(FUSE_LIBRARIES
  NAMES ${fuse_names}
  PATHS /var/services/homes/gbalakrishna/fuse3/fuse-library/lib /usr/lib/fuse-library/lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FUSE DEFAULT_MSG
  FUSE_INCLUDE_DIRS FUSE_LIBRARIES)

mark_as_advanced(
  FUSE_INCLUDE_DIRS FUSE_LIBRARIES)
