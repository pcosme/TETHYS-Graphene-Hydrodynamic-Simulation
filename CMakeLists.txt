cmake_minimum_required(VERSION 3.16)
project(Hydrodynamic_Simulation)

set(CMAKE_CXX_STANDARD 14)

include_directories(src)
include_directories(./src)

find_package(HDF5 REQUIRED COMPONENTS CXX)

find_package(GSL REQUIRED)
link_libraries(GSL::gsl GSL::gslcblas)


set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-unused-parameter -O2")
#set(GSL_FLAGS "-lgsl -lgslcblas -lm")

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -g ")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=leak -g ")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=thread -g ")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=undefined -g ")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=memory -g ")

add_executable(TETHYS_1D
        src/BoundaryLib.cpp
        src/BoundaryLib.h
        src/ElectricLib.cpp
        src/ElectricLib.h
        src/Tethys1DLib.cpp
        src/Tethys1DLib.h
        src/TETHYS_1D_Main_v134.cpp
        src/TethysLib.cpp
        src/TethysLib.h src/TethysMathLib.h src/TethysMathLib.cpp)

add_executable(TETHYS_2D
        src/BoundaryLib.cpp
        src/BoundaryLib.h
        src/ElectricLib.cpp
        src/ElectricLib.h
        src/Tethys2DLib.cpp
        src/Tethys2DLib.h
        src/TETHYS_2D_Main_v221.cpp
        src/TethysLib.cpp
        src/TethysLib.h src/TethysMathLib.h src/TethysMathLib.cpp)

add_executable(TETHYS_ELEC_2D
        src/ElectricLib.cpp
        src/ElectricLib.h
        src/Tethys2DLib.cpp
        src/Tethys2DLib.h
        src/TETHYS_2D_ElectronicAnalysis.cpp
        src/TethysLib.cpp
        src/TethysLib.h src/TethysMathLib.h src/TethysMathLib.cpp)


add_executable(TESTE
        src/BoundaryLib.cpp
        src/BoundaryLib.h
        src/ElectricLib.cpp
        src/ElectricLib.h
        src/Tethys2DLib.cpp
        src/Tethys2DLib.h
        src/DufortFrankel.cpp
        src/TethysLib.cpp
        src/TethysLib.h src/TethysMathLib.h src/TethysMathLib.cpp)


target_include_directories(TETHYS_1D PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(TETHYS_1D PRIVATE ${HDF5_CXX_LIBRARIES})

target_include_directories(TETHYS_2D PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(TETHYS_2D PRIVATE ${HDF5_CXX_LIBRARIES})

target_include_directories(TETHYS_ELEC_2D PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(TETHYS_ELEC_2D PRIVATE ${HDF5_CXX_LIBRARIES})

target_include_directories(TESTE PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(TESTE PRIVATE ${HDF5_CXX_LIBRARIES})