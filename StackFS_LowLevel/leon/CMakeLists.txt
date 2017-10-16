project(StackFS)

cmake_minimum_required(VERSION 2.6)

################################################################################
# The version number.
################################################################################
# The default version number is the latest official build
SET (gatb-tool_VERSION_MAJOR 1)
SET (gatb-tool_VERSION_MINOR 0)
SET (gatb-tool_VERSION_PATCH 0)

# But, it is possible to define another release number during a local build
IF (DEFINED MAJOR)
    SET (gatb-tool_VERSION_MAJOR ${MAJOR})
ENDIF()
IF (DEFINED MINOR)
    SET (gatb-tool_VERSION_MINOR ${MINOR})
ENDIF()
IF (DEFINED PATCH)
    SET (gatb-tool_VERSION_PATCH ${PATCH})
ENDIF()

set (gatb-tool-version ${gatb-tool_VERSION_MAJOR}.${gatb-tool_VERSION_MINOR}.${gatb-tool_VERSION_PATCH})

# However, continuous integration has priority over local compilation
IF (DEFINED JENKINS_TAG)
    SET (gatb-tool-version ${JENKINS_TAG})
ENDIF()

################################################################################
# Define cmake modules directory
################################################################################
SET (GATB_CORE_HOME  ${PROJECT_SOURCE_DIR}/thirdparty/gatb-core/gatb-core)
SET (CMAKE_MODULE_PATH ${GATB_CORE_HOME}/cmake)
SET (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR})

################################################################################
# SUPPORTED KMER SIZES
################################################################################

# One can uncomment this line and set the wanted values
#set (KSIZE_LIST "32   64   96  128  160  192  224  256")

################################################################################
# THIRD PARTIES
################################################################################

# We don't want to install some GATB-CORE artifacts
SET (GATB_CORE_EXCLUDE_TOOLS     1)
SET (GATB_CORE_EXCLUDE_TESTS     1)
SET (GATB_CORE_EXCLUDE_EXAMPLES  1)

# GATB CORE
include (GatbCore)

################################################################################
# TOOLS
################################################################################

find_package (FUSE REQUIRED)
find_package(GLIB2 REQUIRED)
#include_directories (${FUSE_INCLUDE_DIR})
#include_directories (${GLIB_INCLUDE_DIRS})
# We add the compilation options for the library
add_definitions (-D_FILE_OFFSET_BITS=64 -DFUSE_USE_VERSION=30 -g)
add_definitions (${gatb-core-flags})

# We add the gatb-core include directory
include_directories (${gatb-core-includes} ${GLIB2_INCLUDE_DIRS})

# We add the path for extra libraries
link_directories (${gatb-core-extra-libraries-path})

set (PROGRAM_SOURCE_DIR ${PROJECT_SOURCE_DIR}/src)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_FILE_OFFSET_BITS=64")

include_directories (${PROGRAM_SOURCE_DIR})
file (GLOB_RECURSE  ProjectFiles  ${PROGRAM_SOURCE_DIR}/*.cpp)
add_executable(${PROJECT_NAME} ${ProjectFiles})
<<<<<<< HEAD
SET (LIBFUSE_LIB  /usr/local/lib/libfuse3.so)
target_link_libraries(${PROJECT_NAME} ${LIBFUSE_LIB} ${gatb-core-libraries} ${GLIB2_LIBRARIES})
=======
target_link_libraries(${PROJECT_NAME} ${gatb-core-libraries} ${FUSE_LIBRARIES} ${GLIB2_LIBRARIES})
>>>>>>> 88432777b0dd31fa28e88314a24e535fed24f9c5


################################################################################
# DELIVERY
################################################################################

# We set the version number
SET (CPACK_PACKAGE_DESCRIPTION_SUMMARY  "gatb-tool ${PROJECT_NAME}")
SET (CPACK_PACKAGE_VENDOR               "Genscale team (INRIA)")
SET (CPACK_PACKAGE_VERSION_MAJOR        "${gatb-tool_VERSION_MAJOR}")
SET (CPACK_PACKAGE_VERSION_MINOR        "${gatb-tool_VERSION_MINOR}")
SET (CPACK_PACKAGE_VERSION_PATCH        "${gatb-tool_VERSION_PATCH}")
SET (CPACK_PACKAGE_VERSION              "${gatb-tool-version}")  

# We set the kind of archive
SET (CPACK_GENERATOR                    "TGZ")
SET (CPACK_SOURCE_GENERATOR             "TGZ")

# We ignore unwated files for the source archive
SET (CPACK_SOURCE_IGNORE_FILES          
    "^${PROJECT_SOURCE_DIR}/\\.git/"     ;
    "^${PROJECT_SOURCE_DIR}/\\.gitmodules"  ;
    "^${PROJECT_SOURCE_DIR}/\\.gitignore";
    "^${PROJECT_SOURCE_DIR}/build/"  ;
    "^${GATB_CORE_HOME}/\\.cproject" ;
    "^${GATB_CORE_HOME}/\\.git/"     ;
    "^${GATB_CORE_HOME}/\\.project"  ;
    "^${GATB_CORE_HOME}/\\.gitignore";
    "^${GATB_CORE_HOME}/doc/design"  ;
    "^${GATB_CORE_HOME}/DELIVERY.md"
)

# We copy some files to the binary archive
INSTALL (TARGETS ${PROJECT_NAME} DESTINATION bin)
INSTALL (DIRECTORY "${PROJECT_SOURCE_DIR}/data/" DESTINATION "data")
INSTALL (FILES LICENSE README scripts/simple_test.sh DESTINATION bin/..)

# We include the "bin" tag into binary archive file name
set (CPACK_PACKAGE_FILE_NAME  ${PROJECT_NAME}-${CPACK_PACKAGE_VERSION}-bin-${CMAKE_SYSTEM_NAME})

include (CPack)
