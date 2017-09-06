#include <iostream>

#include "matrixsimulator.h"
#include <chrono>
#include "argagg.hpp"

using namespace std::chrono;
using argagg::parser_results;
using argagg::parser;
using std::cerr;
using std::cout;
using std::endl;
using std::ofstream;
using std::ostream;
using std::ostringstream;
using std::string;

int main(int argc, char *argv[])
{

    MatrixSimulator simulator;

    parser argparser {{

                              {
                                      "help", {"-h", "--help"},
                                      "Print help and exit", 0
                              },

                              {

                                      "spectrograph", {"-s", "--spectrograph"},
                                      "path to spectrograph model", 1

                              },

                              {
                                      "blackbody", {"-b", "--blackbody"},
                                      "OPTIONAL: Simulate a blackbody with effective temperature K and magnitude M", 0
                              },

                              {

                                      "temperature", {"-T", "--temperature"},
                                      "OPTIONAL: Input temperature for source", 1

                              },

                              {

                                      "magnitude", {"-m", "--magnitude"},
                                      "OPTIONAL: Input magnitude for source", 1

                              },

                      }};

    // Define our usage text. (Really poor quality)
    ostringstream usage;
    usage
            << argv[0] << " 1.0" << endl
            << endl
            << "Usage: " << argv[0] << " [OPTIONS]... [FILES]..." << endl
            << endl;

    // Use our argument parser to... parse the command line arguments. If there
    // are any problems then just spit out the usage and help text and exit.
    argagg::parser_results args;

    try {
        args = argparser.parse(argc, argv);
    } catch (const std::exception& e) {
        argagg::fmt_ostream fmt(cerr);
        fmt << usage.str() << argparser << endl
            << "Encountered exception while parsing arguments: " << e.what()
            << endl;
        return EXIT_FAILURE;
    }

    // If the help flag was specified then spit out the usage and help text and
    // exit.
    if (args["help"]) {
        argagg::fmt_ostream fmt(cerr);
        fmt << usage.str() << argparser;
        return EXIT_SUCCESS;
    }

    double temp;
    double mag;

    Source * cs = new Source();

    if (args["temperature"]) {

        temp = args["temperature"].as<double>();

    }
    else{
        temp = 6500;
    }

    if (args["magnitude"]) {

        cout<<"Stae";
        mag = args["magnitude"].as<double>();

    }
    else{
        mag = 0;
    }

    if (args["blackbody"]) {

        cout<<"Stae";
        cs = new Blackbody(temp,mag);

    }
    else{

        cs = new Constant(10E-8);

    }



    for (int i=1; i<2; ++i){
        high_resolution_clock::time_point t1 = high_resolution_clock::now();
        simulator.load_spectrograph_model(args["spectrograph"], i, i>1);
        std::cout<< "Fiber " << i << std::endl;
        Telescope Gemini = Telescope();
        GratingEfficiency ge = GratingEfficiency(0.8, simulator.get_blaze(), simulator.get_blaze(), simulator.get_gpmm());
//        EtalonEfficiency ee = EtalonEfficiency(10.,1.,0., 0.95);
//        ConstantEfficiency ge = ConstantEfficiency(1.);
        simulator.add_efficiency(&ge);
        simulator.add_telescope(&Gemini);
//        simulator.add_efficiency(&ee);

//        Blackbody cs = Blackbody(9602., 0.);

//          Constant cs = Constant(10E-8);

//          PhoenixSpectrum cs = PhoenixSpectrum("/data/CppLibs/7125_0_lte03200-5.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_2.fits",
//          "/data/CppLibs/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits",0.45, 0.85, 0.);
//        IdealEtalon cs = IdealEtalon(5., 1., 0., 0.95);
//        cs.set_doppler_shift(-100.);
//        LineList cs = LineList("/home/stuermer/Repos/cpp/EchelleSimulator/laser.txt");
        simulator.add_source(cs);

        simulator.set_wavelength(10000);
//        simulator.set_wavelength(cs.get_wavelength());
//        add_vector_parallel();
//        simulator.prepare_interpolation(1000000);

        //simulator.photon_order(1000000);
          simulator.photon_order(1);
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