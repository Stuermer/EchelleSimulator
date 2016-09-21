//#include <opencv2/imgproc.hpp>
#include <iostream>

#include "matrixsimulator.h"
#include "Slit.h"
#include "efficiency.h"
#include <chrono>
#include "helper.h"
#include "noise.h"
#include "source.h"
#include "PSF.h"
#include "CCD.h"

using namespace std::chrono;

int main(int argc, char *argv[])
{
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    PhoenixSpectrum ps = PhoenixSpectrum("/home/julian/Dissertation/CRIRES-POP/template/lte03400-4.00-1.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits","/home/julian/Dissertation/CRIRES-POP/template/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits", 0.45, 0.7);

    PSF_ZEMAX psfs = PSF_ZEMAX(argv[1]);
    Slit s = Slit(50., 150., 10);
    CCD ccd = CCD(4096, 4096, 3, s.slit_image.type());

    MatrixSimulator simulator;

    simulator.read_transformations(argv[2]);
    simulator.set_wavelength(10000);
    simulator.set_ccd(&ccd);
    simulator.set_slit(&s);
    simulator.set_psfs(&psfs);

    GratingEfficiency ge = GratingEfficiency(0.8, 76., 76., 31.6);
    simulator.add_efficiency(&ge);

    IdealEtalon cs = IdealEtalon(10., 1., 0., 0.9);
    simulator.add_source(&cs);

    simulator.simulate_spectrum();

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>( t2 - t1 ).count()/1000000.;
    std::cout << "Duration: "  << duration << std::endl;

    simulator.save_to_file("../image2.hdf", true, true);

    return 0;
}