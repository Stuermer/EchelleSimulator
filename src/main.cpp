#include <iostream>

#include "matrixsimulator.h"
#include <chrono>

using namespace std::chrono;

int main(int argc, char *argv[])
{

    MatrixSimulator simulator;
//    create_fits_file("../simulations/mdwarf_13bil_0m.fit");

    for (int i=1; i<2; ++i){
        high_resolution_clock::time_point t1 = high_resolution_clock::now();
        simulator.load_spectrograph_model(argv[1], i, i>1);
        std::cout<< "Fiber " << i << std::endl;

        GratingEfficiency ge = GratingEfficiency(0.8, simulator.get_blaze(), simulator.get_blaze(), simulator.get_gpmm());
//        EtalonEfficiency ee = EtalonEfficiency(5.,1.,0., 0.95);
//        ConstantEfficiency ge = ConstantEfficiency(1.);
        simulator.add_efficiency(&ge);
//        simulator.add_efficiency(&ee);

//        Blackbody cs = Blackbody(3500.);
        Constant cs = Constant(1000000.);
//        PhoenixSpectrum cs = PhoenixSpectrum("/data/work/template/7125_0_lte03200-5.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_2.fits",
//        "/data/work/template/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits",0.45, 0.85);
//        IdealEtalon cs = IdealEtalon(5., 1., 0., 0.95);
//        cs.set_doppler_shift(-100.);
//        LineList cs = LineList("/home/stuermer/Repos/cpp/EchelleSimulator/laser.txt");
        simulator.add_source(&cs);

        simulator.set_wavelength(10000);
//        simulator.set_wavelength(cs.get_wavelength());
//        add_vector_parallel();
//        simulator.prepare_interpolation(1000000);

        //simulator.photon_order(1000000);
          simulator.photon_order(1, 1);
//        simulator.photon_order_artifical(50000,0.0002);
//        simulator.simulate_spectrum(false);

        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>( t2 - t1 ).count()/1000000.;
        std::cout << "Total Duration: "  << duration << std::endl;

//        simulator.save_to_fits("../simulations/flat"+std::to_string(i)+".fit", true, false, false);
//        simulator.save_1d_to_fits("../simulations/etalon_noblaze_tri.fit");
    }
    simulator.save_to_fits("../simulations/test.fit", false, false, true);

    return 0;
}