#include <iostream>
#include <chrono>
#include <CCfits/FITS.h>
#include <CCfits/PHDU.h>
#include "argagg.hpp"
#include "helper.h"
#include "matrixsimulator.h"

#define FMT_HEADER_ONLY

#include <fmt/format.h>

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
                                     "OPTIONAL: .csv file(s) with wavelength dependent efficiency values for correct signal scaling. File format is wl[microns];efficiency[fractional]"
                                     "multiple efficiency files can be specified (efficiencies will be multiplied), separate them with ','"
                                     "Make sure, efficiencies cover the full wavelength range of the spectrograph", 1

                             },
                             {

                                     "noblaze", {"--noblaze"},
                                     "OPTIONAL: if set, the simulation will not scale the efficiency with the grating alpha function that was calculated analytically from the spectrograph optical parameters.", 0

                             },
                             {
                                     "blaze", {"--blaze"},
                                     "OPTIONAL: alpha angle in degree. ", 1
                             },
                             {

                                     "etalon", {"--etalon"},
                                     "OPTIONAL: ideal fabry-perot etalon light source"
                                     "general usage: --etalon <d>,<n>,<theta>,<R>,<I>"
                                     "where <d> is the mirror distance in mm, <n> the refractive index between the mirrors"
                                     "<theta> the angle of incidence, <R> the reflectivity of the mirrors, and <I> the constant"
                                     "flux density in [microwatts]/[micrometer] of the light source"
                                     "example usage: --etalon 10,1.,0.,0.94,0.001", 1

                             },
                     }};


    std::ostringstream usage;
    usage
            << argv[0] << " 1.0" << std::endl
            << std::endl
            << "Usage: " << argv[0] << " [OPTIONS]... [FILES]..." << std::endl
            << std::endl;

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

    auto diam = args["telescope"].as<double>(1.);
    auto f_telescope = args["telescope"].as<double>(1.);
    Telescope telescope = Telescope(diam, f_telescope);

    MatrixSimulator simulator(spectrograph, fiber, false);

    auto *cs = new Source();

    if (args["blackbody"]) {
        auto v = args["blackbody"].as<std::string>();
        std::vector<std::string> vv = split_to_vector(v, ',');
        std::cout << fmt::format("Simulating a blockbody with T={:} and V-band magnitude {:}", stod(vv[0]), stod(vv[1]))
                  << std::endl;

        cs = new Blackbody(stod(vv[0]), stod(vv[1]), telescope.get_area());
        source = "blackbody";
    } else if (args["phoenix"]) {
        auto v = args["phoenix"].as<std::string>();
        std::vector<std::string> vv = split_to_vector(v, ',');
        std::cout << "Simulating phoenix spectra with magnitude = " << stod(vv[4]) << std::endl;
        if (download_phoenix(std::stoi(vv[0]), std::stod(vv[3]), std::stod(vv[1]), std::stod(vv[2]),
                             "../data/phoenix_spectra/spectrum.fits") == 0) {

            const std::string &w_file = "../data/phoenix_spectra/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits";

            if (!check_for_file(w_file)) {
                download_wave_grid("../data/phoenix_spectra/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits");
            }
            cs = new PhoenixSpectrum("../data/phoenix_spectra/spectrum.fits",
                                     "../data/phoenix_spectra/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits", stod(vv[4]),
                                     telescope.get_area());

        } else {
            argagg::fmt_ostream fmt(std::cerr);
            fmt << usage.str() << argparser;
            return EXIT_FAILURE;
        }
        source = "phoenix";
    } else if (args["coehlo"]) {

        auto v = args["coehlo"].as<std::string>();
        std::vector<std::string> vv = split_to_vector(v, ',');
        std::cout << "Simulating coehlo spectra with magnitude " << stod(vv[1]) << std::endl;

        cs = new CoehloSpectrum(vv[0], stod(vv[1]), telescope.get_area());
        source = "choehlo";
    } else if (args["custom1"]) {

        auto v = args["custom1"].as<std::string>();
        std::vector<std::string> vv = split_to_vector(v, ',');
        std::cout << "Simulating custom spectra with magnitude " << stod(vv[3]) << std::endl;

        cs = new CustomSpectrum(stod(vv[3]), telescope.get_area(), vv[1], vv[2]);

    } else if (args["custom2"]) {

        auto v = args["custom2"].as<std::string>();
        std::vector<std::string> vv = split_to_vector(v, ',');
        std::cout << "Simulating custom spectra with magnitude = " << stod(vv[2]) << std::endl;

        cs = new CustomSpectrum(std::stod(vv[3]), telescope.get_area(), vv[0], vv[1]);

    } else if (args["linelist"]) {

        auto v = args["linelist"].as<std::string>();
        std::vector<std::string> vv = split_to_vector(v, ',');
        std::cout << "Simulating line list spectra " << vv[0] << std::endl;
        cs = new LineList(vv[0]);

        source = "line list";

    } else if (args["constant"]) {
        auto c_val = args["constant"].as<double>(0.01);
        std::cout << "Simulating constant source with spectral density = " << c_val
                  << " [micro watt] / [micro meter] " << std::endl;

        cs = new Constant(c_val);
        source = "constant";
    } else if (args["etalon"]) {
        auto v = args["etalon"].as<std::string>();
        std::vector<std::string> vv = split_to_vector(v, ',');
        std::cout << fmt::format("Simulating ideal etalon with distance d={:}", stod(vv[0])) << std::endl;
        cs = new IdealEtalon(stod(vv[0]), stod(vv[1]), stod(vv[2]), stod(vv[3]), stod(vv[4]));
        source = "IdealEtalon";
    } else {

        std::cout << "Simulating constant source with spectral density = 0.01 [micro watt] / [micro meter]"
                  << std::endl;
        cs = new Constant(0.01);
        source = "constant";
    }

    auto rv = args["radial_velocity"].as<double>(0.);

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    auto *global_eff = new Efficiency();
    if (args["noblaze"]) {
        global_eff = new ConstantEfficiency(1.);
    } else {
        auto blaze = args["blaze"].as<double>(simulator.get_alpha());
        global_eff = new GratingEfficiency(1., simulator.get_alpha(), blaze, simulator.get_gpmm());
    }
    simulator.add_efficiency(global_eff);

    auto custom_eff = args["efficiency"].as<std::string>("");

    std::vector<Efficiency *> efficienies;

    if (!(custom_eff.empty())) {
        auto v = args["efficiency"].as<std::string>();
        std::vector<std::string> eff_csv_files = split_to_vector(v, ',');
        for (auto ef: eff_csv_files) {
            std::cout << "Loading efficiency curve from " << custom_eff << std::endl;
            efficienies.push_back(new CSVEfficiency(ef));
            simulator.add_efficiency(efficienies.back());
        }
    }

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