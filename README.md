![EchelleSimulator](https://github.com/Stuermer/EchelleSimulator/blob/master/doc/logo.png "Echelle Simulator")


Échelle++, a fast generic échelle simulator.
[![Build Status](https://travis-ci.org/Stuermer/EchelleSimulator.svg?branch=master)](https://travis-ci.org/Stuermer/EchelleSimulator)

The basic idea is that any (fiber-fed) echelle spectrograph can be modelled with a set of wavelength-dependent transformation matrices and point spread functions which describe the spectrographs optics:

First, wavelength-dependent **affine transformation matrices** are extracted from the ZEMAX model of the spectrograph. As the underlying geometric transformations (scaling, rotation, shearing, translation) vary smoothly across an echelle order, these matrices can be interpolated for any intermediate wavelength.

Second, a wavelength-dependent **point spread functions (PSFs)** is applied on the transformed slit images to properly account for optical aberrations. Again, the PSF is only slowly varying across an echelle order, allowing for interpolation at intermediate wavelength.

#### Concept:
---
![Echelle simulation](https://github.com/Stuermer/EchelleSimulator/blob/master/doc/intro.png "Echelle simulation")

---

**Both, the matrices and the PSFs have to be extracted from ZEMAX only once. It is therefore possible to simulate spectra without access to ZEMAX**

---

#### Output:

---
![Echelle spectrum](https://github.com/Stuermer/EchelleSimulator/blob/master/doc/flat.png "Simulated Echelle")

---
Part of a simulated flat spectrum for the MAROON-X spectrograph.

#### Features:

---
 * parallelized C++ code for fast simulations
 * arbitrary 1D spectra
 * arbitrary PSFs
 * arbitrary efficiency models can be applied
 * works with any spectrograph (needs access to ZEMAX model only once)

---

## Prerequisites
 * GCC > 4.9 (or equivalent MSVC), capable of handling C++11 syntax
 * [HDF 5.0 library](https://www.hdfgroup.org/hdf5/)
 * [OpenCV](http://opencv.org/)
 * [Curl](https://curl.haxx.se/libcurl/)
 
## Example usage
```bash
./echellesimulator --spectrograph MaroonX --phoenix 3500,-1.,0.,5.5,1 -r 100 -o mdwarf.fits
```
simulates a phoenix M-dwarf spectrum with the given stellar parameters, and a RV shift of 100m/s.

See
```bash
./echellesimulator -h
``` 
for all available programm arguments.
Check also examples folder for python scripting.

## Documentation
The documentation can be found [here](https://stuermer.github.io/EchelleSimulator).
The documentation is automatically produced by **doxygen**, using [this](https://github.com/Velron/doxygen-bootstrapped) template.

## Contribution
Contributions are welcome! You can help by 
* report bugs
* provide spectrograph models
* make suggestions 
* improve documentation
* improve code
* implement new features
