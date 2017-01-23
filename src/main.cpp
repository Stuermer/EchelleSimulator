#include <iostream>

#include "matrixsimulator.h"
#include <chrono>

using namespace std::chrono;

int main(int argc, char *argv[])
{
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    MatrixSimulator simulator;
    simulator.load_spectrograph_model(argv[1], 1);

    GratingEfficiency ge = GratingEfficiency(0.8, simulator.get_blaze(), simulator.get_blaze(), simulator.get_gpmm());
    simulator.add_efficiency(&ge);

    IdealEtalon cs = IdealEtalon(5., 1., 0., 0.9);
    simulator.add_source(&cs);

    simulator.set_wavelength(10000);
    simulator.simulate_spectrum();

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>( t2 - t1 ).count()/1000000.;
    std::cout << "Duration: "  << duration << std::endl;

    simulator.save_to_fits("../MaroonX.fit", true, false, true);

    return 0;
}