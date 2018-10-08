#include <iostream>
#include <chrono>
#include <CCfits/FITS.h>
#include <CCfits/PHDU.h>
#include "argagg.hpp"
#include "helper.h"
#include "matrixsimulator.h"

using argagg::parser_results;
using argagg::parser;

/**
 * Main entry point for Echelle++
 * See echellesimulator -h for all arguments
 * @param argc program arguments
 * @param argv program arguments
 * @return
 */
int main(int argc, char *argv[]) {

    parser argparser{{

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
                                     "if 1 it adds simulated spectrum to .fits file rather than overwrites it. "
                                     "This can be used for the simulation of multiple fibers/slits.", 1

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
                                     "Example usage: --custom spectra_file,min_w,max_w,magnitude", 1

                             },

                             {

                                     "custom2", {"--custom2"},
                                     "OPTIONAL: Simulate a spectra given by the user."
                                     "Example usage: --custom spectra_file,wave_file,magnitude", 1

                             },

                             {

                                     "linelist", {"--linelist"}, "OPTIONAL: simulates line list", 1

                             },

                             {
                                     "phoenix", {"-p", "--phoenix"},
                                     "OPTIONAL: Simulate a mdwarf phoenix spectra with effective temperature T, magnitude M, log g, metalicity, alpha."
                                     "Check http://phoenix.astro.physik.uni-goettingen.de/?page_id=15 for parameter ranges."
                                     "general usage: --phoenix <T>,<Z>,<alpha>,<log g>,<mag>"
                                     "example usage: --phoenix 3200,-1.,0.,5.5,1", 1
                             },

                             {

                                     "constant", {"-c", "--constant"},
                                     "OPTIONAL: Simulate a constant spectra with density in units [micro watt] / ([micro meter] * [meter]^2) "
                                     "in a wavelength range min_w to max_w [micro meter] "
                                     "general usage:  constant <micro watt>"
                                     "(default: --constant 0.01)", 1

                             },
                             {
                                 "telescope", {"--telescope"},
                                 "OPTIONAL: uses a telescope with diameter <diam> [meter] and a focal length of <fl> [meter] for calculating the photon flux"
                                 "generale usage: telescope <diam>, <fl>"
                                 "example usage: telescope 8.1,128.12 ", 1
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

                             {

                                     "seed", {"--seed"},
                                     "OPTIONAL: random seed used for simulations. If 0, random seed is generated by std::random_device. (default: 0) ", 1

                             },

                             {

                                     "readnoise", {"--readnoise"},
                                     "OPTIONAL: std deviation of readnoise (default:0) ", 1

                             },

                             {

                                     "bias", {"--bias"},
                                     "OPTIONAL: bias level count. (default: 0)", 1

                             },

                             {

                                     "efficiency", {"--efficiency"},
                                     "OPTIONAL: .csv file with wavelength dependent efficiency values for correct signal scaling. File format is wl;efficiency", 1

                             },
                     }};


    // Define our usage text. (Really poor quality)
    std::ostringstream usage;
    usage
            << argv[0] << " 1.0" << std::endl
            << std::endl
            << "Usage: " << argv[0] << " [OPTIONS]... [FILES]..." << std::endl
            << std::endl;

    // Use our argument parser to... parse the command line arguments. If there
    // are any problems then just spit out the usage and help text and exit.
    argagg::parser_results args;

    try {
        args = argparser.parse(argc, argv);
    } catch (const std::exception &e) {
        argagg::fmt_ostream fmt(std::cerr);
        fmt << usage.str() << argparser << std::endl
            << "Encountered exception while parsing arguments: " << e.what()
            << std::endl;
        return EXIT_FAILURE;
    }

    // If the help flag was specified then spit out the usage and help text and
    // exit.
    if (args["help"]) {
        argagg::fmt_ostream fmt(std::cerr);
        fmt << usage.str() << argparser;
        return EXIT_SUCCESS;
    }

    std::string source;
    auto keep = args["keep"].as<bool>(false);
    auto fiber = args["fiber"].as<int>(1);

    auto spectrograph = args["spectrograph"].as<std::string>("MaroonX");
    spectrograph = "../data/spectrographs/" + spectrograph + ".hdf";

    MatrixSimulator simulator(spectrograph, fiber, false);

    auto *cs = new Source();

    if (args["blackbody"]) {
        auto v = args["blackbody"].as<std::string>();
        std::vector<std::string> vv = split_to_vector(v, ',');
        std::cout << "Simulating a blackbody with T = " << stod(vv[0]) << " and magnitude K = " << stod(vv[0]) << std::endl;

        cs = new Blackbody(stod(vv[0]), stod(vv[1]));
        source = "blackbody";
    } else if (args["phoenix"]) {
        auto v = args["phoenix"].as<std::string>();
        std::vector<std::string> vv = split_to_vector(v, ',');
        std::cout << "Simulating phoenix spectra with magnitude = " << stod(vv[4]) << std::endl;
        if (download_phoenix(std::stoi(vv[0]), std::stod(vv[3]), std::stod(vv[1]), std::stod(vv[2]), "../data/phoenix_spectra/spectrum.fits") == 0) {

            const std::string &w_file = "../data/phoenix_spectra/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits";

            if (!check_for_file(w_file)) {
                download_wave_grid("../data/phoenix_spectra/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits");
            }
            cs = new PhoenixSpectrum("../data/phoenix_spectra/spectrum.fits",
                                     "../data/phoenix_spectra/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits", stod(vv[4]));

        } else {
            argagg::fmt_ostream fmt(std::cerr);
            fmt << usage.str() << argparser;
            return EXIT_FAILURE;
        }
        source = "phoenix";
    } else if (args["coehlo"]) {

        auto v = args["choehlo"].as<std::string>();
        std::vector<std::string> vv = split_to_vector(v, ',');
        std::cout << "Simulating coehlo spectra with magnitude = " << stod(vv[1]) << std::endl;

        cs = new CoehloSpectrum(vv[0], stod(vv[1]));
        source = "choehlo";
    } else if (args["custom1"]) {

        auto v = args["custom1"].as<std::string>();
        std::vector<std::string> vv = split_to_vector(v, ',');
        std::cout << "Simulating coehlo spectra with magnitude = " << stod(vv[3]) << std::endl;

        cs = new CustomSpectrum(vv[0], stod(vv[1]), stod(vv[2]), stod(vv[3]));

    } else if (args["custom2"]) {

        auto v = args["custom2"].as<std::string>();
        std::vector<std::string> vv = split_to_vector(v, ',');
        std::cout << "Simulating coehlo spectra with magnitude = " << stod(vv[2]) << std::endl;

        cs = new CustomSpectrum(vv[0], vv[1], stod(vv[2]));

    } else if (args["linelist"]) {

        auto v = args["linelist"].as<std::string>();
        std::vector<std::string> vv = split_to_vector(v, ',');
        std::cout << "Simulating line list spectra with scaling factor = " << stod(vv[1]) << std::endl;
        cs = new LineList(vv[0], stod(vv[1]));

        source = "line list";

    } else if (args["constant"]) {
        auto v = args["constant"].as<std::string>();
        std::vector<std::string> vv = split_to_vector(v, ',');
        std::cout << "Simulating constant source with spectral density = " << stod(vv[0])
             << " [micro watt] / ([micro meter] * [meter]^2)" << std::endl;

        cs = new Constant(stod(vv[0]), simulator.get_minimum_wavelength(), simulator.get_maximum_wavelength());
        source = "constant";
    } else {

        std::cout << "Simulating constant source with spectral density = 1 [micro watt] / ([micro meter] * [meter]^2)"
             << std::endl;
        cs = new Constant();
        source = "constant";
    }

    auto rv = args["radial_velocity"].as<double>(0.);

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    auto diam = args["telescope"].as<double>(1.);
    auto f_telescope = args["telescope"].as<double>(1.);
    Telescope telescope = Telescope(diam, f_telescope);
    GratingEfficiency ge = GratingEfficiency(0.8, simulator.get_blaze(), simulator.get_blaze(), simulator.get_gpmm());
    auto ef = args["efficiency"].as<std::string>("");

    auto *eff = new Efficiency();;
    if (!(ef.empty())) {
        std::cout << "Loading efficiency curve from " << ef << std::endl;
        eff = new CSVEfficiency(ef);
        simulator.add_efficiency(eff);
    }

    simulator.add_efficiency(&ge);
    simulator.set_telescope(&telescope);

    cs->set_doppler_shift(rv);
    simulator.set_source(cs);

    // in case of 'normal' continuous spectrum
    if (!cs->is_list_like()) {
        simulator.set_wavelength(10000);
    } else {
        simulator.set_wavelength(cs->get_wavelength());
        std::cout << "Running LineList Test" << std::endl;
    }

    auto t = args["integration_time"].as<double>(1.);
    auto seed = args["seed"].as<unsigned long>(0);
    auto readnoise = args["readnoise"].as<double>(0);
    auto bias = args["bias"].as<double>(0);

    simulator.simulate(t, seed);
    simulator.add_background(bias, readnoise, seed);

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000000.;
    std::cout << "Total Duration: " << duration << std::endl;

    auto path = args["output"].as<std::string>("test.fit");
    if (path.find('/') == std::string::npos) {
        simulator.save_to_fits("../simulations/" + path, false, !keep);
        std::string filename = "../simulations/" + path;
        std::unique_ptr<CCfits::FITS> pFits;
        pFits.reset(new CCfits::FITS(filename, CCfits::Write));

        try {
            pFits->pHDU().addKey("EXPTIME", t, "exposure time");
            pFits->pHDU().addKey("Spectrograph", spectrograph, "Used spectrograph model");
        }
        catch (...) {

            std::cout << "keys already exists - skipping";
        };
        pFits->pHDU().addKey("Fiber_" + std::to_string(fiber), source, "Source for Fiber " + std::to_string(fiber));

    } else {
        simulator.save_to_fits(path, false, !keep);
    }

    return 0;
}