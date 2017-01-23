project(echellesimulator)


cmake_minimum_required(VERSION 2.8)
find_package( OpenCV REQUIRED )
#add_definitions(-DUSE_GPU)
set(CUDA_USE_STATIC_CUDA_RUNTIME OFF)

set (LINK_LIBS ${LINK_LIBS} ${OpenCV_LIBS})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/src/CMakeModules/")

#INCLUDE("${CMAKE_SOURCE_DIR}/src/CMakeModules/FindCCFITS.cmake")
#set (LINK_LIBS ${LINK_LIBS} ${CCFITS_LIBRARIES})

find_package(CCFITS REQUIRED)
# Install stage
INSTALL_FILES("${CMAKE_ROOT}/Modules" FILES "${CMAKE_SOURCE_DIR}/FindCCOIFITS.cmake")

FIND_LIBRARY(CCFITS_LIB CCfits /usr/lib/x86_64-linux-gnu)

find_package(HDF5 COMPONENTS HL REQUIRED)

include_directories(SYSTEM ${HDF5_INCLUDE_DIRS} ${HDF5_HL_INCLUDE_DIR})
include_directories(SYSTEM $OPENCV_INCLUDE_DIRS)
include_directories(src)
include_directories(include)
set (LINK_LIBS ${LINK_LIBS} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} ${CCFITS_LIB} ${CFITSIO_LIBRARY})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
 set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_C_FLAGS} -fopenmp -O3")
#set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_C_FLAGS} -O3")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")

add_executable(echellesimulator src/main.cpp src/helper.cpp src/source.cpp src/efficiency.cpp src/Slit.cpp src/matrixsimulator.cpp src/noise.cpp src/PSF.cpp src/CCD.cpp src/hdf5opencv.cpp include/csv_reader.h)
target_link_libraries(echellesimulator ${LINK_LIBS})
