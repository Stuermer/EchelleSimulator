#include <iostream>

#include "matrixsimulator.h"
#include <chrono>
#include <CCfits/FITS.h>
#include <CCfits/PHDU.h>
#include "argagg.hpp"
#include "helper.h"

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

                              "coehlo", {"--coehlo"},
                              "OPTIONAL: Simulate a solar coehlo spectra from a file."
                              "Check http://specmodels.iag.usp.br/fits_search/?refer=s_coelho05 for available files."
                              "Example usage: --coehlo file_path,magnitude", 1

                      },

                      {

                              "custom1", {"--custom1"},
                              "OPTIONAL: Simulate a spectra given by the user."
                              "Example usage: --custom spectra_file,min_w,max_w,magnitude",1

                      },

                      {

                              "custom2", {"--custom2"},
                              "OPTIONAL: Simulate a spectra given by the user."
                              "Example usage: --custom spectra_file,wave_file,magnitude",1

                      },

                      {

                              "linelist", {"--linelist"},"Hello",1

                      },

                      {
                              "phoenix", {"-p", "--phoenix"},
                              "OPTIONAL: Simulate a mdwarf phoenix spectra with effective temperature T, magnitude M, log g, metalicity alpha."
                                  "Check http://phoenix.astro.physik.uni-goettingen.de/?page_id=15 for parameter ranges - intermediate values will be rounded to available spectra."
                              "example usage: --phoenix 3200,1.,-5.5,0.,1", 1
                              // The explicit command for our current setup m=0
                              //--spectrograph MaroonX -p 7000,1.50,-0.5,1.0,0.0
                      },

                      {

                              "constant", {"-c", "--constant"},
                              "Optional: Simulate a constant spectra with density in units [micro watt] / ([micro meter] * [meter]^2) "
                              "in a wavelength range min_w to max_w [micro meter] (default: --constant 1 0 1)", 1

                      },

                      {

                              "radial_velocity", {"-r", "--radial-velocity"},
                              "OPTIONAL: radial velocity shift [m/s] (default: 0) ", 1

                      },
                      {

                              "integration_time", {"-t", "--integration-time"},
                              "OPTIONAL: integration time of the spectrograph [s] (default: 1) ", 1

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

    std::string source;
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
        cout<<"Simulating a blackbody with T = "<< stod(vv[0])  <<" and magnitude K = " << stod(vv[0]) << endl;

        cs = new Blackbody(stod(vv[0]) , stod(vv[1]));
        source = "blackbody";
    }
    else if (args["phoenix"]) {
        auto v = args["phoenix"].as<string>();
        std::vector<std::string> vv = split(v, ',');
        cout<<"Simulating phoenix spectra with magnitude = "<< stod(vv[4]) << endl;

        if(download_phoenix(vv[0], vv[1], vv[2], vv[3]) == 0){

            const std::string& w_file = "../data/phoenix_spectra/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits";

            if(!check_for(w_file)) {
                download_wave_grid("../data/phoenix_spectra/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits");
            }
            //cs = new CustomSpectrum("../data/phoenix_spectra/test.fits", "../data/phoenix_spectra/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits", stod(vv[4]));
//            cout<<"Downloaded";
            cs = new PhoenixSpectrum("../data/phoenix_spectra/test.fits", "../data/phoenix_spectra/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits", stod(vv[4]));

        }
        else{

            argagg::fmt_ostream fmt(cerr);
            fmt << usage.str() << argparser;
            return EXIT_FAILURE;

        }
        source = "phoenix";

    }
    else if (args["coehlo"]) {

        auto v = args["choehlo"].as<string>();
        std::vector<std::string> vv = split(v, ',');
        cout<<"Simulating coehlo spectra with magnitude = "<< stod(vv[1]) << endl;

        cs = new CoehloSpectrum(vv[0],stod(vv[1]));

    }
    else if (args["custom1"]) {

        auto v = args["custom1"].as<string>();
        std::vector<std::string> vv = split(v, ',');
        cout<<"Simulating coehlo spectra with magnitude = "<< stod(vv[3]) << endl;

        cs = new CustomSpectrum(vv[0],stod(vv[1]),stod(vv[2]),stod(vv[3]));

    }
    else if (args["custom2"]) {

        auto v = args["custom2"].as<string>();
        std::vector<std::string> vv = split(v, ',');
        cout<<"Simulating coehlo spectra with magnitude = "<< stod(vv[2]) << endl;

        cs = new CustomSpectrum(vv[0],vv[1],stod(vv[2]));

    }
    else if (args["linelist"]){

        auto v = args["linelist"].as<string>();
        std::vector<std::string> vv = split(v, ',');
        cout<<"Simulating line list spectra with scaling factor = "<< stod(vv[1]) << endl;
        cs = new LineList(vv[0], stod(vv[1]));

        simulator.mode = false;

        source = "line list";

    }
    else if (args["constant"]) {
        auto v = args["constant"].as<string>();
        std::vector<std::string> vv = split(v, ',');
        cout<<"Simulating constant source with spectral density = "<< stod(vv[0]) << " [micro watt] / ([micro meter] * [meter]^2)" << endl;

        cs = new Constant(stod(vv[0]),stod(vv[1]),stod(vv[2]));
        source = "constant";
    }
    else{

        cout<<"Simulating constant source with spectral density = 1 [micro watt] / ([micro meter] * [meter]^2)" << endl;
        cs = new Constant();
        source = "constant";
    }

    auto rv = args["radial_velocity"].as<double>(0.);

    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    Telescope Gemini = Telescope();
    GratingEfficiency ge = GratingEfficiency(0.8, simulator.get_blaze(), simulator.get_blaze(), simulator.get_gpmm());
//    ConstantEfficiency ge = ConstantEfficiency(1.);

    simulator.add_efficiency(&ge);
    simulator.set_telescope(&Gemini);

    cs->set_doppler_shift(rv);
    simulator.set_source(cs);

//    simulator.set_wavelength(10000);
    if(cs -> mode == 1){

        simulator.set_wavelength(10000);

    }
    else {

        simulator.set_wavelength(cs -> get_wavelength());

    }

    auto t = args["integration_time"].as<double>(1.);
    simulator.simulate(t);

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>( t2 - t1 ).count()/1000000.;
    std::cout << "Total Duration: "  << duration << std::endl;

    auto path = args["output"].as<std::string>("test.fit");
    if (path.find("/") == std::string::npos) {
        simulator.save_to_fits("../simulations/" + path, false, false, true);
        std::string filename = "../simulations/" + path;
        std::vector<std::string> * keys;
        std::auto_ptr<CCfits::FITS> pFits(0);
        pFits.reset( new CCfits::FITS(filename, CCfits::Write));

        try {
            pFits->pHDU().addKey("EXPTIME", t, "exposure time");
            pFits->pHDU().addKey("Spectrograph", spectrograph, "Used spectrograph model");
        }
        catch(...) {

                std::cout<< "keys already exists - skipping";
        };
        pFits->pHDU().addKey("Fiber_"+std::to_string(fiber), source, "Source for Fiber "+std::to_string(fiber));

    }
    else {
        simulator.save_to_fits(path);
    }

    return 0;
}