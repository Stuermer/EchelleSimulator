project(echellesimulator)

add_executable( echellesimulator cvmatandqimage.cpp efficiency.cpp slit.cpp matrixsimulator.cpp main.cpp )
cmake_minimum_required(VERSION 2.8)
find_package( OpenCV REQUIRED )
find_package(CUDA QUIET REQUIRED)
install(TARGETS echellesimulator RUNTIME DESTINATION bin)

set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
# set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
target_link_libraries( echellesimulator ${OpenCV_LIBS} )


# Find the QtWidgets library
find_package(Qt5Widgets)
find_package(Qt5Core)
# Use the Widgets module from Qt 5.
target_link_libraries(echellesimulator Qt5::Widgets)

