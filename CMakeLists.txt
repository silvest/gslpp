cmake_minimum_required(VERSION 3.5)

project(gslpp)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

set(CMAKE_CXX_FLAGS "-Wall ${GSL_CFLAGS} ${ROOT_CFLAGS}")

include_directories(${BOOST_INC})

file(GLOB srcs "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
add_library(gslpp ${LIBTYPE} ${srcs})

file(GLOB headers "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h")
INSTALL(FILES ${headers} DESTINATION include/HEPfit COMPONENT header)
