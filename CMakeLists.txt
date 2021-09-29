project(echellesimulator)
cmake_minimum_required(VERSION 3.0)
set(CMAKE_CXX_STANDARD 11)
include(ExternalProject)

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

# Download fmt
ExternalProject_Add(fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt.git
        GIT_TAG 3e75ad9822980e41bc591938f26548f24eb88907
        INSTALL_COMMAND ""
        UPDATE_COMMAND ""
        BUILD_COMMAND ""
        )
ExternalProject_Get_Property(fmt source_dir)
set(FMT_INCLUDE_DIR ${source_dir}/include)
INCLUDE_DIRECTORIES(${FMT_INCLUDE_DIR})

# Download argagg
ExternalProject_Add(argagg
        GIT_REPOSITORY https://github.com/vietjtnguyen/argagg.git
        GIT_TAG 4c8c86180cfafb1448f583ed0973da8c2f559dd6
        INSTALL_COMMAND ""
        UPDATE_COMMAND ""
        BUILD_COMMAND ""
        )
ExternalProject_Get_Property(argagg source_dir)
MESSAGE(${source_dir})
set(ARGAGG_INCLUDE_DIR ${source_dir}/include/argagg)
MESSAGE(${ARGAGG_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${ARGAGG_INCLUDE_DIR})


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")

file(GLOB SOURCES
        include/*.h
        src/*.cpp
        )

list(FILTER SOURCES EXCLUDE REGEX ".*main.cpp$")

add_library(echellelib STATIC ${SOURCES})
set (LINK_LIBS ${LINK_LIBS} ${HDF5_LIBRARIES} ${CCFITS_LIB} ${CFITSIO} ${HDF5_HL_LIBRARIES} curl)
MESSAGE(${LINK_LIBS})
target_link_libraries(echellelib ${LINK_LIBS})


add_executable(echellesimulator src/main.cpp)
add_dependencies(echellesimulator fmt)
add_dependencies(echellesimulator argagg)
target_include_directories(echellesimulator PRIVATE ${FMT_INCLUDE_DIR})
target_link_libraries(echellesimulator ${LINK_LIBS} echellelib)

if (${TESTS})
 add_subdirectory(tests)
endif ()
