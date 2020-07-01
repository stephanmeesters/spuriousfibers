cmake_minimum_required(VERSION 2.8)
 
if(POLICY CMP0020)
  cmake_policy(SET CMP0020 NEW)
endif()

PROJECT(SpuriousFibers)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

# Check for C++14 support
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#find_package(OpenMP)
#if (OPENMP_FOUND)
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#	add_compile_definitions(WITH_OPENMP)
#endif()






find_library(OpenMP_LIBRARY
    NAMES omp
)

find_path(OpenMP_INCLUDE_DIR
    omp.h
)

mark_as_advanced(OpenMP_LIBRARY OpenMP_INCLUDE_DIR)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OpenMP DEFAULT_MSG 
    OpenMP_LIBRARY OpenMP_INCLUDE_DIR)

if (OpenMP_FOUND)
    set(OpenMP_LIBRARIES ${OpenMP_LIBRARY})
    set(OpenMP_INCLUDE_DIRS ${OpenMP_INCLUDE_DIR})
    set(OpenMP_COMPILE_OPTIONS -Xpreprocessor -fopenmp)

    add_library(OpenMP::OpenMP SHARED IMPORTED)
    set_target_properties(OpenMP::OpenMP PROPERTIES
        IMPORTED_LOCATION ${OpenMP_LIBRARIES}
        INTERFACE_INCLUDE_DIRECTORIES "${OpenMP_INCLUDE_DIRS}"
        INTERFACE_COMPILE_OPTIONS "${OpenMP_COMPILE_OPTIONS}"
    )
endif()










include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ./lib)

file(GLOB CXX_FILES *.cxx)

add_executable(SpuriousFibers ${CXX_FILES})

target_link_libraries(SpuriousFibers ${VTK_LIBRARIES} ${OpenMP_LIBRARY})



























