#ifndef MATRIXSIMULATOR_H
#define MATRIXSIMULATOR_H

#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <memory>

#include "helper.h"
#include "spline.h"
#include "efficiency.h"
#include "source.h"
#include "PSF.h"
#include "CCD.h"
#include "Slit.h"
#include "telescope.h"
#include <array>

struct point2d{
    double x,y;
};

struct raw_transformation {
    int order;
    double wavelength;
    std::array<float,6> transformation_matrix;
    std::vector<double> decomposed_matrix;

};

struct spectrograph_information {
    double blaze;
    double gpmm;
};

typedef struct raw_transformation raw_transformation;

class MatrixSimulator {
public:
    /**
     * Load spectrograph model from HDF file
     * @param path path to HDF file containing spectrograph model
     * @param fiber_number fiber to select
     * @param keep_ccd if true it assumes that a CCD has already been added to the spectrograph model. It keeps it and
     * adds the new simulations to it. This can be used to simplify the simulation of multiple fibers.
     */
    MatrixSimulator(std::string path, int fiber_number, bool keep_ccd);

    /**
     * Set wavelength grid on which the input spectrum will be interpolated.
     * @param N number of wavelength per order
     */
    void set_wavelength(int N);

    void set_wavelength(std::vector<double> wavelength);

    /**
     * Adds an efficiency profile to the simulator.
     * @param eff
     */
    void add_efficiency(Efficiency *eff);

    /**
     * Sets telescope of the spectrograph.
     * @param telescope
     */
    void set_telescope(Telescope *telescope);

    /**
     * Sets the spectral source of the current simulation.
     * @param src
     */
    void set_source(Source *src);

    /**
     * Save simulated echelle image to an HDF file.
     * @param filename filename
     * @param bleed bleed overexposed pixel TODO: not implemented yet
     * @param overwrite True to overwrite existing file TODO: not working yet
     */
    void save_to_hdf(std::string filename, bool downsample = true, bool bleed = true, bool overwrite = false);

    /**
     *
     * @param filename filename
     * @param bleed bleed overexposed pixel TODO: not implemented yet
     * @param overwrite True to overwrite existing file TODO: not working yet
     */
    void save_to_fits(std::string filename, bool downsample = true, bool bleed = true, bool overwrite = false);

    /**
     * Save 1d spectra in fits file.
     * This function saves the 1 dimensional spectrum which was used for the simulation in the fits file.
     * It is therefore a perfectly reduced spectrum.
     * @param filename path to the fits file
     */
    void save_1d_to_fits(std::string filename);

    /**
     * Returns blaze angle in degrees
     * @return blaze angle
     */
    double get_blaze();

    /**
     * Returns grove lines per mm
     * @return grove lines per mm
     */
    double get_gpmm();

    /**
     * Returns fiber number of current simulation.
     * @return current fiber number
     */
    int get_fiber_number();

    /**
     * Simulate echelle spectrum.
     * @param t integration time
     * @return
     */
    int simulate(double t);

    /**
     * Returns the minimum wavelength supported by the spectrograph [microns]
     * @return minimum wavelength [microns]
     */
    double get_minimum_wavelength();

    /**
    * Returns the maximum wavelength supported by the spectrograph [microns]
    * @return maximum wavelength [microns]
    */
    double get_maximum_wavelength();

    bool mode = true;

private:
    /**
     * Load spectrograph model from HDF file
     * @param path path to HDF file containing spectrograph model
     * @param fiber_number fiber to select
     * @param keep_ccd if true it assumes that a CCD has already been added to the spectrograph model. It keeps it and
     * adds the new simulations to it. This can be used to simplify the simulation of multiple fibers.
     */
    void load_spectrograph_model(std::string path, int fiber_number, bool keep_ccd = false);

    /**
     *
     */
    void calc_splines();
    /**
 * Get affine transformation matrix at specific wavelength and order
 * @param order echelle diffraction order
 * @param wavelength wavelength [micron]
 * @return 2x3 affine transformation matrix
 */
    std::array<float,6> get_transformation_matrix(int order, double wavelength);

    /**
     * Get affine transformation matrix, but use lookup tables for speedup.
     * In good approximation the parameters shear, rotation, scale_x and scale_y will not vary quickly.
     * Use lookup tables for them for speedup.
     * @param o echelle diffraction order
     * @param wavelength wavelength [micron]
     * @return 2x3 affine transformation matrix
     */
    inline std::array<float,6> get_transformation_matrix_lookup(int o, double wavelength);

    void set_efficiencies(std::vector<Efficiency *> &efficiencies);

    void set_ccd(CCD *ccd);

    void set_slit(Slit *slit);

    void set_psfs(PSF *psfs);

    void prepare_sources(std::vector<Source *> sources);
    void prepare_psfs(int N);
    void prepare_matrix_lookup(int N);

    std::vector<int> orders;
    cv::Mat img;
    int fiber_number;
    int n_orders;
    int min_order;
    int max_order;
    double wavelength_limit_max = 0.; // will be overwritten by load_spectrograph model
    double wavelength_limit_min = 100.; // will be overwritten by load_spectrograph model

    std::map<int, std::vector<raw_transformation> > raw_transformations;

    std::vector<Efficiency *> efficiencies;
    std::vector<Source *> sources;
    Telescope telescope;
    CCD *ccd;
    PSF *psfs;
    Slit *slit;

    std::vector<std::vector<double>> sim_wavelength;
    std::vector<std::vector<std::array<float,6> >> sim_matrices;
    std::vector<std::vector<double>> sim_efficiencies;
    std::vector<std::vector<float>> sim_spectra;
    std::vector<std::vector<float>> sim_spectra_time_efficieny;

    std::vector<float> sim_total_efficiency_per_order;
    std::vector<std::vector<cv::Mat>> sim_psfs;
    std::vector<std::vector<double>> sim_psfs_wavelength;
    std::vector<double> sim_psfs_dwavelength;

    std::vector<std::vector<double>> sim_p;
    std::vector<std::vector<double>> sim_q;
    std::vector<std::vector<double>> sim_r;
    std::vector<std::vector<double>> sim_phi;
    std::vector<std::vector<double>> sim_tx;
    std::vector<std::vector<double>> sim_ty;

    std::vector<std::vector<double>> sim_m00;
    std::vector<std::vector<double>> sim_m01;
    std::vector<std::vector<double>> sim_m10;
    std::vector<std::vector<double>> sim_m11;

    std::vector<std::vector<double>> sim_matrix_wavelength;
    std::vector<double> sim_matrix_dwavelength;

    std::vector<tk::spline> tr_p;
    std::vector<tk::spline> tr_r;
    std::vector<tk::spline> tr_q;
    std::vector<tk::spline> tr_phi;
    std::vector<tk::spline> tr_tx;
    std::vector<tk::spline> tr_ty;

    spectrograph_information spec_info;
};

#endif // MATRIXSIMULATOR_H
