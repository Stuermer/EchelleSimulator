Installation {#installation}
============

## Prerequisites
 * GCC>4.6, cabable of handling C++11 syntax (MSVC is untested, but should work)
 * [HDF 5.0 library](https://www.hdfgroup.org/hdf5/)
 * [OpenCV 2.4](http://opencv.org/)
 * [CCFits](http://heasarc.gsfc.nasa.gov/fitsio/ccfits/)
 
When compiling openCV, make sure you have activated the appropriate compiler flags.
For full execution speed, use release compilation mode. If you want to use CUDA,
make sure you have activated it.
### Installation on ubuntu 16.04 (or similar):
Install hdf library:
    
    sudo apt-get install libhdf5-dev

Install CCFits:

    sudo apt-get install libccfits-dev

Install OpenCV:
    
    wget https://github.com/Itseez/opencv/archive/2.4.13.zip
    unzip 2.4.13.zip
    cd opencv-2.4.13
    mkdir release
    cd release
    cmake ../
    # you can also use cmake-gui for an easy configuration of the compiler flags
    make
    sudo make install

Install Echelle++:
    
    git clone https://github.com/Stuermer/EchelleSimulator.git
    cd EchelleSimulator
    cmake .
    make
## 