cmake_minimum_required(VERSION 2.8)
 
if(POLICY CMP0020)
  cmake_policy(SET CMP0020 NEW)
endif()

PROJECT(SpuriousFibers)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ./lib)

file(GLOB CXX_FILES *.cxx)

add_executable(SpuriousFibers ${CXX_FILES})

target_link_libraries(SpuriousFibers ${VTK_LIBRARIES})



























