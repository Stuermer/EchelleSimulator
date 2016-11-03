![EchelleSimulator](https://github.com/Stuermer/EchelleSimulator/blob/master/doc/logo.png "Echelle Simulator")


Échelle++, a fast generic échelle simulator.


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
 * CUDA support (not fully functional yet)
 * arbitrary 1D spectra
 * arbitrary PSFs
 * arbitrary efficiency modells can be applied
 * works with any spectrograph (needs access to ZEMAX modell only once)

---

## Prerequisites
 * GCC > 4.6 (or equivalent MSVC), capable of handling C++11 syntax
 * [HDF 5.0 library](https://www.hdfgroup.org/hdf5/)
 * [OpenCV 2.4](http://opencv.org/)
 
## Example usage
```c++
#include <iostream>
#include "matrixsimulator.h"

int main(int argc, char *argv[])
{
PSF_ZEMAX psfs = PSF_ZEMAX(argv[1]); // Load ZEMAX simulated psfs from file
Slit s = Slit(50., 150., 10);        // define slit size
CCD ccd = CCD(4096, 4096, 3, s.slit_image.type()); //define detector

MatrixSimulator simulator;  //create instance of simulator

simulator.read_transformations(argv[2]);  //read in ZEMAX simulated transformation matrices
simulator.set_wavelength(10000); // set wavelength steps per order
simulator.set_ccd(&ccd); // set simulator ccd
simulator.set_slit(&s); //set simulator slit
simulator.set_psfs(&psfs); //set simulator psfs

GratingEfficiency ge = GratingEfficiency(0.8, 76., 76., 31.6); //Echelle Grating efficiency
simulator.add_efficiency(&ge); //add efficiency profile to simulator. More profiles can be added

IdealEtalon cs = IdealEtalon(10., 1., 0., 0.9); // spectral source ideal etalon
simulator.add_source(&cs); //add source to simulator

simulator.simulate_spectrum(); //simulate echelle spectra 

simulator.save_to_file("../image2.hdf", true, true); //save echelle spectrum to file

return 0;
}
```

generates a 2D Fabry-Perot etalon spectrum.

On a mordern PC this takes about 3s - 5s.

## Documentation
The documentation can be found [here](https://stuermer.github.io/EchelleSimulator).
The documentation is automatically produced by **doxygen**, using [this](https://github.com/Velron/doxygen-bootstrapped) template.

## Contribution
Contributions are welcome! You can help by 
* report bugs
* provide spectroraph models
* make suggestions 
* improve documentation
* improve code
* implement new features
