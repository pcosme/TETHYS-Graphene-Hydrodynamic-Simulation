cmake_minimum_required(VERSION 3.16)
project(Hydrodynamic_Simulation)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib")

enable_testing()

include_directories(src/)
include_directories(src/includes)

find_package(HDF5 REQUIRED COMPONENTS CXX)
include_directories(${HDF5_INCLUDE_DIRS})

#find_package(GSL REQUIRED)
#link_libraries(GSL::gsl GSL::gslcblas)
#set(GSL_FLAGS "-lgsl -lm")
#set(BLAS_FLAGS "-llapack -lblas -llapacke -lm")

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-unused-parameter -O2")
FIND_PACKAGE(OpenMP REQUIRED)
if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_C_FLAGS}")
endif()


#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -g ")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=leak -g ")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=thread -g ")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=undefined -g ")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=memory -g ")

add_executable(TETHYS_1D src/TETHYS_1D_Main_v134.cpp)
add_executable(TETHYS_2D src/TETHYS_2D_Main_v250.cpp)
add_executable(TETHYS_ELEC_2D src/TETHYS_2D_ElectronicAnalysis.cpp)


add_test(test1D TETHYS_1D 20 10 0.003 0.002 0.001 0.03  0.0 1 1)
add_test(test2D TETHYS_2D 20 10 0.003 0.002 0.001 0.03  0.0 1 .5)
add_test(testElec2D TETHYS_ELEC_2D hdf5_2D_S=20.00vF=10.00vis=0.00odd=0.001l=0.00wc=0.00.h5')

file(GLOB SRC_LIB_FILES src/*Lib.cpp)
add_library(tethys_library STATIC ${SRC_LIB_FILES})

find_library (TETHYS_LIB tethys_library lib)
target_link_libraries(TETHYS_1D LINK_PUBLIC ${TETHYS_LIB})
target_link_libraries(TETHYS_2D LINK_PUBLIC ${TETHYS_LIB})
target_link_libraries(TETHYS_ELEC_2D LINK_PUBLIC ${TETHYS_LIB})

target_include_directories(TETHYS_1D PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(TETHYS_1D PRIVATE ${HDF5_CXX_LIBRARIES})

target_include_directories(TETHYS_2D PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(TETHYS_2D PRIVATE ${HDF5_CXX_LIBRARIES})

target_include_directories(TETHYS_ELEC_2D PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(TETHYS_ELEC_2D PRIVATE ${HDF5_CXX_LIBRARIES})
