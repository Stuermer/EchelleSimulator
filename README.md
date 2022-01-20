![EchelleSimulator](https://github.com/Stuermer/EchelleSimulator/blob/master/docs/plots/logo.png "Echelle Simulator")


Echelle++, a fast generic spectrum simulator.
[![Build Status](https://travis-ci.org/Stuermer/EchelleSimulator.svg?branch=master)](https://travis-ci.org/Stuermer/EchelleSimulator)
[![Doc Status](https://readthedocs.org/projects/echellesimulator/badge/?version=latest)](https://echellesimulator.readthedocs.io/en/latest/)
[![DOI](https://zenodo.org/badge/68249504.svg)](https://zenodo.org/badge/latestdoi/68249504)



Echelle++ is a simulation tool, to generate realistic 2D spectra, in particular cross-dispersed echelle spectra.
It allows to simulate arbitrary spectra for any fiber-fed or slit spectrograph, where a model file
is available. Optical aberrations are treated accurately, the simulated spectra include photon and read-out noise.  

:stop_sign: --------------------------------------------------- :stop_sign: ---------------------------------------------------  :stop_sign:
## ECHELLE++ is no longer actively developed. Please check out [PyEchelle](https://gitlab.com/Stuermer/pyechelle)  instead.
Even though I might implement bug fixes, I moved the project over to python (with more features such as CUDA support) for easier deployment. 

:stop_sign: --------------------------------------------------- :stop_sign: ---------------------------------------------------  :stop_sign:

### Example usage
```bash
./echellesimulator --spectrograph MaroonX --phoenix 3500,-1.,0.,5.5,1 -r 100 -o mdwarf.fits
```
simulates a phoenix M-dwarf spectrum with the given stellar parameters, and a RV shift of 100m/s for the MAROON-X spectrograph.
The output looks similar to:

![Echelle spectrum](https://github.com/Stuermer/EchelleSimulator/blob/master/docs/plots/mdwarf.jpg "Simulated Echelle")

See
```bash
./echellesimulator -h
``` 
for all available program arguments.
Check also the examples folder for python scripting.


### Features:

 * parallel C++ code for fast simulations
 * arbitrary 1D input spectra
 * arbitrary PSFs
 * arbitrary efficiency models can be applied
 * works with any spectrograph (needs access to ZEMAX model only once)
 * currently provided spectrographs: MAROON-X, NEID, VeloceRosso


### Concept:
The basic idea is that any (fiber-fed) echelle spectrograph can be modelled with a set of wavelength-dependent 
transformation matrices and point spread functions which describe the spectrographs optics:

First, wavelength-dependent **affine transformation matrices** are extracted from the ZEMAX model of the spectrograph. 
As the underlying geometric transformations (scaling, rotation, shearing, translation) vary smoothly across an echelle 
order, these matrices can be interpolated for any intermediate wavelength.

Second, a wavelength-dependent **point spread functions (PSFs)** is applied on the transformed slit images to properly 
account for optical aberrations. Again, the PSF is only slowly varying across an echelle order, allowing for 
interpolation at intermediate wavelength.

![Echelle simulation](https://github.com/Stuermer/EchelleSimulator/blob/master/docs/plots/intro.png "Echelle simulation")

**Both, the matrices and the PSFs have to be extracted from ZEMAX only once. It is therefore possible to simulate 
spectra without access to ZEMAX**

For more information see [here](https://echellesimulator.readthedocs.io/en/latest/).

### How to use Echelle++
There are two ways of using Echelle++: Building from source, or using the [docker](https://www.docker.com/) image.
#### Building
##### Prerequisites
 * GCC > 4.9 (or equivalent MSVC), capable of handling C++11 syntax
 * CMake >=3.0
 * [CCFits](https://heasarc.gsfc.nasa.gov/fitsio/ccfits/)
 * [CFITSIO](https://heasarc.gsfc.nasa.gov/fitsio/fitsio.html)
 * [HDF 5.0 library](https://www.hdfgroup.org/hdf5/)
 * [Curl](https://curl.haxx.se/libcurl/)
 * [fmt](https://github.com/fmtlib/fmt) Header-only library will be downloaded automatically by CMake

Install the required 3rd party packages. Make sure they are compiled with the same compiler version to avoid runtime issues.
Run cmake and make to build Echelle++.

For Linux, see [here](https://stuermer.github.io/EchelleSimulator/installation.html) for a full installation guide.

#### Docker
For convenience, a docker image is provided that runs on any platform as long as docker is set up correctly (see [here](https://www.docker.com/get-started)).
After docker is installed a simple
```bash
docker run -v /path/to/output_directory:/home/simulations stuermer/echellesimulator
``` 
will download the latest version of Echelle++ and run it with the given arguments.
So, a 
```bash
docker run -v /path/to/output_directory:/home/simulations stuermer/echellesimulator --spectrograph MaroonX --phoenix 3500,-1.,0.,5.5,1 -r 100 -o mdwarf.fits
```
 will start the simulation as shown above and save the output *mdwarf.fit* on your local folder */path/to/output_directory*
 
See [here](https://echellesimulator.readthedocs.io/en/latest/installation.html) for platform dependent considerations.
 
### Documentation
The package documentation can be found [here](https://echellesimulator.readthedocs.io/en/latest/index.html).

### Contribution
Contributions are welcome! You can help by 
* report bugs
* provide spectrograph models
* make suggestions 
* improve documentation
* improve code
* implement new features
