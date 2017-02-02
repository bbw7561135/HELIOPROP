cmake_minimum_required(VERSION 2.6)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(pugixml STATIC
    src/pugixml.cpp
)

SET_TARGET_PROPERTIES( pugixml PROPERTIES COMPILE_FLAGS -DWITH_PUGIXML)

