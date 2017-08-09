project(echellesimulator)


cmake_minimum_required(VERSION 2.8)
#include_directories(/home/stuermer/)
find_package( OpenCV REQUIRED HINTS /data/CppLibs/opencv/)
#add_definitions(-DUSE_GPU)
set(CUDA_USE_STATIC_CUDA_RUNTIME OFF)

set(LINK_LIBS ${LINK_LIBS} ${OpenCV_LIBS})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/src/CMakeModules/")

#INCLUDE("${CMAKE_SOURCE_DIR}/src/CMakeModules/FindCCFITS.cmake")
#set (LINK_LIBS ${LINK_LIBS} ${CCFITS_LIBRARIES})

#find_package(CCFITS REQUIRED HINTS /home/stuermer/libs/CCfits/usr/local/lib)
#Install stage
INSTALL_FILES("${CMAKE_ROOT}/Modules" FILES "${CMAKE_SOURCE_DIR}/FindCCOIFITS.cmake")

FIND_LIBRARY(CCFITS_LIB CCfits /data/CppLibs/CCFits/usr/local/lib)
include_directories(/data/CppLibs/CCFits/usr/local/include)
include_directories(/data/CppLibs/opencv/include/opencv)

#find_package(HDF5 COMPONENTS HL REQUIRED)
find_package(HDF5 COMPONENTS C CXX HL REQUIRED)

FIND_LIBRARY(CFITSIO cfitsio /usr/local/lib64)
# set(HDF5_INCLUDE_DIRS /home/stuermer/Repos/hdf5-1.8.18/c++/src)
#set(HDF5_HL_INCLUDE_DIR /home/stuermer/Repos/hdf5-1.8.18/hl/src)

#include_directories(SYSTEM ${HDF5_INCLUDE_DIRS} ${HDF5_HL_INCLUDE_DIR})

include_directories(SYSTEM $OPENCV_INCLUDE_DIRS)
include_directories(src)
include_directories(include)
#include_directories(/home/stuermer/libs/include)
include_directories(/data/CppLibs/eigen-eigen-67e894c6cd8f)
#FIND_LIBRARY(HDF5_LIBRARIES hdf5 /home/stuermer/Repos/hdf5-1.8.18/src/.libs)
#FIND_LIBRARY(HDF5_HL_LIBRARIES hdf5_hl /home/stuermer/Repos/hdf5-1.8.18/hl/src/.libs)
#link_directories( ${HDF5_LIBRARY_DIRS} )
include_directories( ${HDF5_INCLUDE_DIRS} )

set (LINK_LIBS ${LINK_LIBS} ${HDF5_LIBRARIES} ${CCFITS_LIB} ${CFITSIO} ${HDF5_HL_LIBRARIES})
MESSAGE(${LINK_LIBS})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -O3")
#set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -fopenmp -O3")
#set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_C_FLAGS} -O3")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")

add_executable(echellesimulator src/main.cpp src/helper.cpp src/source.cpp src/efficiency.cpp src/Slit.cpp src/matrixsimulator.cpp src/noise.cpp src/PSF.cpp src/CCD.cpp src/hdf5opencv.cpp include/csv_reader.h src/telescope.cpp include/telescope.h)
target_link_libraries(echellesimulator ${LINK_LIBS})
