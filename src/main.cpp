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

    //download_phoenix("ftp://phoenix.astro.physik.uni-goettingen.de/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-1.0.Alpha=+0.80/lte06000-4.50-1.0.Alpha=+0.80.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits",
                    // "../data/phoenix_spectra/test.fits");

    parser argparser {{

                      {
                              "help", {"-h", "--help"},
                              "Print help and exit", 0
                      },

                      {

                              "spectrograph", {"-s", "--spectrograph"},
                              "name of spectrograph (default: MaroonX), name has to match filename", 1

                      },

                      {

                              "fiber", {"-f", "--fiber"},
                              "fiber number to use for simulations. Starts at 1. (default: 1) ", 1

                      },

                      {

                              "keep", {"-k", "--keep-ccd"},
                              "if 1 it assumes that a CCD has already been added to the spectrograph model. It keeps it and "
                               "adds the new simulations to it. This can be used to simplify the simulation of multiple fibers.", 1

                      },

                      {
                              "blackbody", {"-b", "--blackbody"},
                              "OPTIONAL: Simulate a blackbody with effective temperature K and magnitude M. example usage: --blackbody 3500,1.0 ", 1
                      },

                      {
                              "phoenix", {"-p", "--phoenix"},
                              "OPTIONAL: Simulate a mdwarf phoenix spectra with effective temperature T, magnitude M, log g, metalicity alpha."
                                  "Check http://phoenix.astro.physik.uni-goettingen.de/?page_id=15 for parameter ranges - intermediate values will be rounded to available spectra."
                              "example usage: --phoenix 3200,1.,-5.5,0.,1", 1
                              // The explicit command for our current setup m=0
                              //-r 0 --spectrograph ../data/MaroonX.hdf  -p "/data/CppLibs/7125_0_lte03200-5.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_2.fits","/data/CppLibs/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits",0.
                      },

                      {

                              "constant", {"-c", "--constant"},
                              "Optional: Simulate a constant spectra with density in units [micro watt] / ([micro meter] * [meter]^2) "
                              "in a wavelength range min_w to max_w [micro meter] (default: --constant 1 0 1)", 1

                      },

                      {

                              "radial_velocity", {"-r", "--radial-velocity"},
                              "OPTIONAL: radial velocity shift in m/s (default: 0) ", 1

                      },
                      {

                              "integration_time", {"-t", "--integration-time"},
                              "OPTIONAL: integration time of the spectrograph (default: 1) ", 1

                      },
                      {

                              "output", {"-o", "--output"},
                              "OPTIONAL: path of the output fits file. If only a filename is given, the image will be saved in "
                                      "../simulations/filename. Otherwise it's assumed the path is absolute (default: test.fit) ", 1

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

    auto keep = args["keep"].as<double>(0);
    auto fiber = args["fiber"].as<double>(1);

    auto spectrograph = args["spectrograph"].as<std::string>("MaroonX");
    spectrograph = "../data/spectrographs/" + spectrograph + ".hdf";

    MatrixSimulator simulator(spectrograph, fiber, keep);

    double temp;
    double mag;

    Source * cs = new Source();

    if (args["blackbody"]) {
        auto v = args["blackbody"].as<string>();
        std::vector<std::string> vv = split(v, ',');
        cout<<"Simulating a blackbody with T="<< stod(vv[0])  <<" and magnitude K=" << stod(vv[0]) << endl;

        cs = new Blackbody(stod(vv[0]) , stod(vv[1]));
    }
    else if (args["phoenix"]) {
        auto v = args["phoenix"].as<string>();
        std::vector<std::string> vv = split(v, ',');
        cout<<"Simulating phoenix spectra with magnitude="<< stod(vv[4]) << endl;

        download_phoenix(vv[0], vv[1], vv[2], vv[3], "../data/phoenix_spectra/test.fits");

        cs = new PhoenixSpectrum("../data/phoenix_spectra/test.fits", "/data/CppLibs/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits", 0.45, 0.85, stod(vv[4])); //0.45 and 0.85 are hardcoded wavelength range parameters (they need to be coded out)
    }
    else if (args["constant"]) {
        auto v = args["constant"].as<string>();
        std::vector<std::string> vv = split(v, ',');
        cout<<"Simulating constant source with spectral density="<< stod(vv[0]) << " [micro watt] / ([micro meter] * [meter]^2)" << endl;

        cs = new Constant(stod(vv[0]),stod(vv[1]),stod(vv[2]));
    }
    else{

        cout<<"Simulating constant source with spectral density=1 [micro watt] / ([micro meter] * [meter]^2)" << endl;
        cs = new Constant();
    }

    auto rv = args["radial_velocity"].as<double>(0.);

    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    Telescope Gemini = Telescope();
//    GratingEfficiency ge = GratingEfficiency(0.8, simulator.get_blaze(), simulator.get_blaze(), simulator.get_gpmm());
    ConstantEfficiency ge = ConstantEfficiency(1.);

    simulator.add_efficiency(&ge);
    simulator.set_telescope(&Gemini);

    cs->set_doppler_shift(rv);
    simulator.set_source(cs);

    simulator.set_wavelength(10000);

    auto t = args["integration_time"].as<double>(1.);
    simulator.simulate(t);

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>( t2 - t1 ).count()/1000000.;
    std::cout << "Total Duration: "  << duration << std::endl;

    auto path = args["output"].as<std::string>("test.fit");
    if (path.find("/") == std::string::npos)
        simulator.save_to_fits("../simulations/" + path, false, false, true);
    else
        simulator.save_to_fits(path);

    return 0;
}