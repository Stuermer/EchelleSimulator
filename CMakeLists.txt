project(echellesimulator)

cmake_minimum_required(VERSION 3.0)
find_package( OpenCV REQUIRED HINTS /data/CppLibs/opencv/)
set(CUDA_USE_STATIC_CUDA_RUNTIME OFF)
set(LINK_LIBS ${LINK_LIBS} ${OpenCV_LIBS})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/src/CMakeModules/")

INSTALL_FILES("${CMAKE_ROOT}/Modules" FILES "${CMAKE_SOURCE_DIR}/FindCCOIFITS.cmake")

FIND_LIBRARY(CCFITS_LIB CCfits /data/CppLibs/CCFits/usr/local/lib)
include_directories(/data/CppLibs/CCFits/usr/local/include)
include_directories(/data/CppLibs/opencv/include/opencv)

find_package(HDF5 COMPONENTS C CXX HL REQUIRED)

FIND_LIBRARY(CFITSIO cfitsio /usr/local/lib64)

include_directories(SYSTEM $OPENCV_INCLUDE_DIRS)
include_directories(src)
include_directories(include)

# include eigen library
find_package(Eigen3)
if(NOT EIGEN3_FOUND)
    # Fallback to cmake_modules
    find_package(cmake_modules REQUIRED)
    find_package(Eigen REQUIRED)
    set(EIGEN3_INCLUDE_DIRS ${EIGEN_INCLUDE_DIRS})
    set(EIGEN3_LIBRARIES ${EIGEN_LIBRARIES})  # Not strictly necessary as Eigen is head only
    # Possibly map additional variables to the EIGEN3_ prefix.
else()
    set(EIGEN3_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR})
endif()
INCLUDE_DIRECTORIES (${EIGEN3_INCLUDE_DIR})

set (LINK_LIBS ${LINK_LIBS} ${HDF5_LIBRARIES} ${CCFITS_LIB} ${CFITSIO} ${HDF5_HL_LIBRARIES} curl)
MESSAGE(${LINK_LIBS})
#for optimized speed use following flags

IF(CMAKE_BUILD_TYPE MATCHES Debug)
    message("Debug build.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -std=c++11 ")
ELSEIF(CMAKE_BUILD_TYPE MATCHES Release)
    message("Release build.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -O3")
    set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -std=c++11 -fopenmp -O3")
ELSEIF(CMAKE_BUILD_TYPE MATCHES RelNoParallel)
    message("No Parallelism.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")
    set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -std=c++11 -O3")
ELSEIF(CMAKE_BUILD_TYPE MATCHES RelWithDebInfo)
    message("profiling build")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -g")
    set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -std=c++11 -O3 -g")
ENDIF()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")

add_executable(echellesimulator src/main.cpp src/helper.cpp src/source.cpp src/efficiency.cpp src/Slit.cpp src/matrixsimulator.cpp src/noise.cpp src/PSF.cpp src/CCD.cpp src/hdf5opencv.cpp include/csv_reader.h src/telescope.cpp include/telescope.h)
target_link_libraries(echellesimulator ${LINK_LIBS})
