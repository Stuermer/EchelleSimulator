Installation {#installation}
============

## Prerequisites
 * GCC>4.8, cabable of handling C++11 syntax (MSVC is untested, but should work)
 * [HDF 5.0 library](https://www.hdfgroup.org/hdf5/)
 * [OpenCV](http://opencv.org/)
 * [CCFits](http://heasarc.gsfc.nasa.gov/fitsio/ccfits/)
 * [Curl](https://curl.haxx.se/libcurl/)
 
When compiling openCV, make sure you have activated the appropriate compiler flags.
For full execution speed, use release compilation mode. If you want to use CUDA,
make sure you have activated it.
### Installation on ubuntu 16.04 (or similar):
Install required libraries:
    
    sudo apt-get install libhdf5-dev libccfits-dev libopencv-dev libcurl4-openssl-dev

Install Echelle++:
    
    git clone https://github.com/Stuermer/EchelleSimulator.git
    cd EchelleSimulator
    cmake .
    make 