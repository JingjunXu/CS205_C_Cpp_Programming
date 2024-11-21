cmake_minimum_required(VERSION 3.12)

add_definitions(-DWITH_AVX2)
add_definitions(-mavx)
add_definitions(-O3)

set(CMAKE_CXX_STANDARD 11)

project(pro)

ADD_EXECUTABLE(pro main.cpp mul.cpp)

set(BLA_VENDER OpenBLAS)
find_package(BLAS REQUIRED)
if (BLAS_FOUND)
    #message(STATUS "BLAS library found: ${BLAS_LIBRARIES}")
    target_link_libraries(pro /usr/lib/x86_64-linux-gnu/libopenblas.so)
    message("openblas Found")
else()
    message("Can't find Openblas library")
endif()

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message("OpenMP found.")
    target_link_libraries(pro OpenMP::OpenMP_CXX)
else()
    message("Can't find OpenMP")
endif()