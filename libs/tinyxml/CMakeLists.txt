cmake_minimum_required(VERSION 2.6)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(tinyxml STATIC
    src/tinystr.cpp
    src/tinyxml.cpp
    src/tinyxmlerror.cpp
    src/tinyxmlparser.cpp
)

SET_TARGET_PROPERTIES( tinyxml PROPERTIES COMPILE_FLAGS -DTIXML_USE_STL )
