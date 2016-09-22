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

    //PhoenixSpectrum ps = PhoenixSpectrum("/home/julian/Dissertation/CRIRES-POP/template/lte03400-4.00-1.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits","/home/julian/Dissertation/CRIRES-POP/template/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits", 0.45, 0.7);

    PSF_ZEMAX psfs1 = PSF_ZEMAX(argv[1]);
    PSF_ZEMAX psfs2 = PSF_ZEMAX(argv[2]);
    PSF_ZEMAX psfs3 = PSF_ZEMAX(argv[3]);
    PSF_ZEMAX psfs4 = PSF_ZEMAX(argv[4]);
    PSF_ZEMAX psfs5 = PSF_ZEMAX(argv[5]);

    Slit s = Slit(50., 150., 10);
    CCD ccd = CCD(4096, 4096, 3, s.slit_image.type());

    MatrixSimulator simulator;

    simulator.read_transformations(argv[6]);
    simulator.set_wavelength(10000);
    simulator.set_ccd(&ccd);
    simulator.set_slit(&s);
    simulator.set_psfs(&psfs1);

    GratingEfficiency ge = GratingEfficiency(0.8, 76., 76., 31.6);
    simulator.add_efficiency(&ge);

    IdealEtalon cs = IdealEtalon(10., 1., 0., 0.9);
    simulator.add_source(&cs);

    simulator.simulate_spectrum();


    MatrixSimulator simulator2;

    simulator2.read_transformations(argv[7]);
    simulator2.set_wavelength(10000);
    simulator2.set_ccd(&ccd);
    simulator2.set_slit(&s);
    simulator2.set_psfs(&psfs2);
    simulator2.add_efficiency(&ge);
    simulator2.add_source(&cs);
    simulator2.simulate_spectrum();

    MatrixSimulator simulator3;

    simulator3.read_transformations(argv[8]);
    simulator3.set_wavelength(10000);
    simulator3.set_ccd(&ccd);
    simulator3.set_slit(&s);
    simulator3.set_psfs(&psfs3);
    simulator3.add_efficiency(&ge);
    simulator3.add_source(&cs);
    simulator3.simulate_spectrum();


    MatrixSimulator simulator4;

    simulator4.read_transformations(argv[9]);
    simulator4.set_wavelength(10000);
    simulator4.set_ccd(&ccd);
    simulator4.set_slit(&s);
    simulator4.set_psfs(&psfs4);
    simulator4.add_efficiency(&ge);
    simulator4.add_source(&cs);
    simulator4.simulate_spectrum();


    MatrixSimulator simulator5;

    simulator5.read_transformations(argv[10]);
    simulator5.set_wavelength(10000);
    simulator5.set_ccd(&ccd);
    simulator5.set_slit(&s);
    simulator5.set_psfs(&psfs5);
    simulator5.add_efficiency(&ge);
    simulator5.add_source(&cs);
    simulator5.simulate_spectrum();

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>( t2 - t1 ).count()/1000000.;
    std::cout << "Duration: "  << duration << std::endl;

    simulator.save_to_file("../image2.hdf", true, true);

    return 0;
}