project(echellesimulator)

cmake_minimum_required(VERSION 3.0)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/src/CMakeModules/")
INSTALL_FILES("${CMAKE_ROOT}/Modules" FILES "${CMAKE_SOURCE_DIR}/FindCCOIFITS.cmake")

FIND_LIBRARY(CCFITS_LIB CCfits ${CCFITS_HINT})
include_directories(${CCFITS_INCLUDE_HINT})

find_package(HDF5 COMPONENTS C CXX HL REQUIRED)
FIND_LIBRARY(CFITSIO cfitsio ${CFITSIO_HINT})

include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/include)
message(${CMAKE_SOURCE_DIR})

INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIRS})

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

file(GLOB SOURCES
        include/*.h
        src/*.cpp
        )

add_executable(echellesimulator src/main.cpp src/helper.cpp src/source.cpp src/efficiency.cpp src/Slit.cpp src/matrixsimulator.cpp src/noise.cpp src/PSF.cpp src/CCD.cpp include/csv_reader.h src/telescope.cpp ${SOURCES})
target_link_libraries(echellesimulator ${LINK_LIBS})
