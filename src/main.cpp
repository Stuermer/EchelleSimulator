#include <iostream>

#include "matrixsimulator.h"
#include <chrono>

using namespace std::chrono;

int main(int argc, char *argv[])
{
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    MatrixSimulator simulator;
    create_fits_file("../simulations/et4.fit");

    for (int i=1; i<2; ++i){
        simulator.load_spectrograph_model(argv[1], i, i>1);

//        GratingEfficiency ge = GratingEfficiency(0.8, simulator.get_blaze(), simulator.get_blaze(), simulator.get_gpmm());
        ConstantEfficiency ge = ConstantEfficiency(1.);
        simulator.add_efficiency(&ge);

//        Blackbody cs = Blackbody(3500);
//        Constant cs = Constant(1.);
//        PhoenixSpectrum cs = PhoenixSpectrum("/data/work/template/7125_0_lte03200-5.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_2.fits",
//        "/data/work/template/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits",0.45, 0.75);
//        IdealEtalon cs = IdealEtalon(10., 1., 0., 0.95);
//        cs.set_doppler_shift(-100.);
        LineList cs = LineList("/home/stuermer/Repos/cpp/EchelleSimulator/laser.txt");
        simulator.add_source(&cs);

//        simulator.set_wavelength(1000000);
        simulator.set_wavelength(cs.get_wavelength());
        simulator.photon_order(1000000);
//        simulator.simulate_spectrum(false);

        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>( t2 - t1 ).count()/1000000.;
        std::cout << "Duration: "  << duration << std::endl;

//        simulator.save_to_fits("../MaroonX.fit", true, false, true);
//        simulator.save_1d_to_fits("../simulations/etalon_noblaze_tri.fit");
    }
    simulator.save_to_fits("../simulations/et4.fit", false, false, true);

    return 0;
}