# EchelleSimulator

EchelleSimulator is a tool to simulate realistic arbitrary 2D Echelle spectra. 


The basic idea is that any (fiber-fed) echelle spectrograph can be modelled with a set of wavelength-dependent transformation matrices and point spread functions which describe the spectrographs optics:

First, wavelength-dependent **affine transformation matrices** are extracted from the ZEMAX model of the spectrograph. As the underlying geometric transformations (scaling, rotation, shearing, translation) vary smoothly across an echelle order, these matrices can be interpolated for any intermediate wavelength.

Second, a wavelength-dependent **point spread functions (PSFs)** is applied on the transformed slit images to properly account for optical aberrations. Again, the PSF is only slowly varying across an echelle order, allowing for interpolation at intermediate wavelength.

![Echelle simulation](https://github.com/Stuermer/EchelleSimulator/blob/master/doc/intro.png "Echelle simulation")

**Both, the matrices and the PSFs have to be extracted from ZEMAX only once. It is therefore possible to simulate spectra without access to ZEMAX**


## Prerequisites
 * GCC > 4.6 (or equivalent MSVC), capable of handling C++11 syntax
 * [HDF 5.0 library](https://www.hdfgroup.org/hdf5/)
 * [OpenCV 2.4](http://opencv.org/)
 
## Example usage
  todo
 
