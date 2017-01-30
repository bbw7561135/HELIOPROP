project(helioprop CXX)

cmake_minimum_required(VERSION 2.6)

set( CMAKE_VERBOSE_MAKEFILE off )
set( CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" )
set( CMAKE_CXX_COMPILER "g++" )
set( CMAKE_CXX_FLAGS "-O2 -ftree-vectorize" )
set( HELIOPROP_EXTRA_INCLUDES )
set( HELIOPROP_EXTRA_LIBRARIES )

include_directories( "${CMAKE_SOURCE_DIR}/include" )

set( HELIOPROP_SRCS 
     src/bfield.cpp
     src/helioprop.cpp
     src/input.cpp
     src/particles.cpp
)

# TINYXML (provided)
add_subdirectory(libs/tinyxml)
list(APPEND HELIOPROP_EXTRA_INCLUDES libs/tinyxml/include)
list(APPEND HELIOPROP_EXTRA_LIBRARIES tinyxml)

# GSL (external)
find_package(GSL)
if(GSL_FOUND)
  set(HAVE_LIBGSL 1)
  list(APPEND HELIOPROP_EXTRA_INCLUDES ${GSL_INCLUDE_DIR})
  list(APPEND HELIOPROP_EXTRA_LIBRARIES ${GSL_LIBRARIES})
else(GSL_FOUND)
  message(FATAL_ERROR "Require GSL. Set GSL_HOME")
endif(GSL_FOUND)

# ROOT (external)
#find_package(ROOT)
#if(ROOT_FOUND)
#	add_definitions(-DHELIOPROP_HAVE_ROOT)
#	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CFLAGS}")
#	list(APPEND HELIOPROP_EXTRA_LIBRARIES ${ROOT_LIBS})
#	set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${ROOT_LIBS}")
#	set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${ROOT_LIBS}")
#	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ROOT_LIBS}")
#endif(ROOT_FOUND)

# MAKE EXECUTABLE
include_directories (include ${HELIOPROP_EXTRA_INCLUDES})

add_executable(HELIOPROP ${HELIOPROP_SRCS} )
target_link_libraries (HELIOPROP ${HELIOPROP_EXTRA_LIBRARIES})

# Show summary.
message (STATUS "CMAKE_SYSTEM .......... = ${CMAKE_SYSTEM}")
message (STATUS "BUILD_SHARED_LIBS ..... = ${BUILD_SHARED_LIBS}")
message (STATUS "CMAKE_CXX_COMPILER .... = ${CMAKE_CXX_COMPILER}")
message (STATUS "CMAKE_CXX_FLAGS ....... = ${CMAKE_CXX_FLAGS}")
message (STATUS "CMAKE_INCLUDE_PATH .... = ${CMAKE_INCLUDE_PATH}")