Example Usage {#usage}
=============
```
#include <iostream>
#include "matrixsimulator.h"

int main(int argc, char *argv[])
{
MatrixSimulator simulator;  //create instance of simulator

simulator.load_spectrograph_model(argv[1], 1);  //read in ZEMAX simulated transformation matrices
simulator.set_wavelength(10000); // set wavelength steps per order

GratingEfficiency ge = GratingEfficiency(0.8, 76., 76., 31.6); //Echelle Grating efficiency
simulator.add_efficiency(&ge); //add efficiency profile to simulator. More profiles can be added

IdealEtalon cs = IdealEtalon(10., 1., 0., 0.9); // spectral source ideal etalon
simulator.add_source(&cs); //add source to simulator

simulator.simulate_spectrum(); //simulate echelle spectra 

simulator.save_to_fits("../image2.hdf", true, true, true); //save echelle spectrum to file

return 0;
}
```

generates a 2D Fabry-Perot etalon spectrum.

On a mordern PC this should take <3s.