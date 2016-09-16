project(echellesimulator)


cmake_minimum_required(VERSION 2.8)
find_package( OpenCV REQUIRED )
#add_definitions(-DUSE_GPU)
set(CUDA_USE_STATIC_CUDA_RUNTIME OFF)

set (LINK_LIBS ${LINK_LIBS} ${OpenCV_LIBS})
find_package(CUDA QUIET REQUIRED)

find_package(HDF5 COMPONENTS HL REQUIRED)
include_directories(SYSTEM ${HDF5_INCLUDE_DIRS} ${HDF5_HL_INCLUDE_DIR})
include_directories(SYSTEM $OPENCV_INCLUDE_DIRS)
set (LINK_LIBS ${LINK_LIBS} ${HDF5_LIBRARIES})

set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O2")
# set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")

add_executable( echellesimulator helper.cpp source.cpp efficiency.cpp slit.cpp matrixsimulator.cpp main.cpp noise.cpp noise.h PSF.cpp PSF.h)
target_link_libraries(echellesimulator ${LINK_LIBS})
