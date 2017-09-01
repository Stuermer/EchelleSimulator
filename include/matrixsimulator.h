#ifndef MATRIXSIMULATOR_H
#define MATRIXSIMULATOR_H

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <memory>

#include "helper.h"
#include "spline.h"
#include "efficiency.h"
#include "source.h"
#include "PSF.h"
#include "CCD.h"
#include "Slit.h"
#include "telescope.h"
#include <Eigen/Dense>

struct point2d{
    double x,y;
};

struct raw_transformation {
    int order;
    double wavelength;
    Matrix23f transformation_matrix;
    std::vector<double> decomposed_matrix;

};

struct spectrograph_information {
    double blaze;
    double gpmm;
};

typedef struct raw_transformation raw_transformation;

class MatrixSimulator {
public:
    MatrixSimulator();

    /**
     * Get affine transformation matrix at specific wavelength and order
     * @param order echelle diffraction order
     * @param wavelength wavelength in micron
     * @return 2x3 affine transformation matrix
     */
    Matrix23f get_transformation_matrix(int order, double wavelength);

    Matrix23f get_transformation_matrix_lookup(int o, double wavelength);

    void calc_splines();

    void set_wavelength(int N);

    void set_wavelength(std::vector<double> wavelength);

    void calc_sim_matrices();

#ifdef USE_GPU
    cv::gpu::GpuMat transform_slit(cv::gpu::GpuMat& slit_image, cv::Mat& transformation_matrix, double weight);
    int simulate_order(int order, cv::gpu::GpuMat& slit_image, cv::gpu::GpuMat& output_image, bool aberrations);
    void simulate_spectrum(cv::gpu::GpuMat& slit_image);
#endif

    cv::Mat transform_slit(cv::Mat &slit_image, cv::Mat &transformation_matrix, double weight);

    int simulate_order(int order, cv::Mat &slit_image, cv::Mat &output_image, bool aberrations);

    void simulate_spectrum(bool aberrations);

    void set_efficiencies(std::vector<Efficiency *> &efficiencies);

    void add_efficiency(Efficiency *eff);

    void add_telescope(Telescope *telescope);

    void add_source(Source *src);

    /**
     * Load spectrograph model from HDF file
     * @param path path to HDF file containing spectrograph model
     * @param fiber_number fiber to select
     * @param keep_ccd if true it assumes that a CCD has already been added to the spectrograph model. It keeps it and
     * adds the new simulations to it. This can be used to simplify the simulation of multiple fibers.
     */
    void load_spectrograph_model(std::string path, int fiber_number, bool keep_ccd = false);

    void set_order_range(int min_order, int max_order);

    std::vector<int> orders;

    void set_ccd(CCD *ccd);

    void set_slit(Slit *slit);

    void set_psfs(PSF *psfs);

    int raw_n;

    void prepare_sources(std::vector<Source *> sources);

    void save_to_hdf(std::string filename, bool downsample = true, bool bleed = true, bool overwrite = false);

    void save_to_fits(std::string filename, bool downsample = true, bool bleed = true, bool overwrite = false);

    /**
     * Save 1d spectra in fits file.
     * This function saves the 1 dimensional spectrum which was used for the simulation in the fits file.
     * It is therefore a perfectly reduced spectrum.
     * @param filename path to the fits file
     */
    void save_1d_to_fits(std::string filename);

    void transformation_to_file(std::string filename);

    CCD *ccd;

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

    int get_fiber_number();

    int photon_order(double t);

    int photon_order_artifical(int N_photons, double dl);

    void prepare_psfs(int N);

    void prepare_matrix_lookup(int N);

private:
    cv::Mat img;
    int fiber_number;
    int n_orders;
    int min_order;
    int max_order;
    double wavelength_limit_max = 0.; // will be overwritten by load_spectrograph model
    double wavelength_limit_min = 100.; // will be overwritten by load_spectrograph model

    std::map<int, std::vector<raw_transformation> > raw_transformations;
//    std::vector<std::vector<raw_transformation>> raw_transformations;

    std::vector<Efficiency *> efficiencies;
    std::vector<Source *> sources;
    Telescope telescope;

//    std::map<int, std::vector<double> > sim_wavelength;
    std::vector<std::vector<double>> sim_wavelength;
//    std::map<int, std::vector<Matrix23f> > sim_matrices;
    std::vector<std::vector<Matrix23f>> sim_matrices;
//    std::map<int, std::vector<double> > sim_efficiencies;
    std::vector<std::vector<double>> sim_efficiencies;
//    std::map<int, std::vector<float> > sim_spectra;
    std::vector<std::vector<float>> sim_spectra;
//    std::map<int, std::vector<float> > sim_spectra_time_efficieny;
    std::vector<std::vector<float>> sim_spectra_time_efficieny;

//    std::map<int, float > sim_total_efficiency_per_order;
    std::vector<float> sim_total_efficiency_per_order;
//    std::map<int, std::vector<cv::Mat> > sim_psfs;
    std::vector<std::vector<cv::Mat>> sim_psfs;
//    std::map<int, std::vector<double> > sim_psfs_wavelength;
    std::vector<std::vector<double>> sim_psfs_wavelength;
//    std::map<int, double > sim_psfs_dwavelength;
    std::vector<double> sim_psfs_dwavelength;

    std::vector<std::vector<double>> sim_p;
    std::vector<std::vector<double>> sim_q;
    std::vector<std::vector<double>> sim_r;
    std::vector<std::vector<double>> sim_phi;
    std::vector<std::vector<double>> sim_matrix_wavelength;
    std::vector<double> sim_matrix_dwavelength;



    //    std::map<int, std::map<float, std::vector<cv::Mat> > psf_lookup;
//    std::vector<double> IP_wl;
//    std::vector<float> INTERPOLATION_eff;
//    std::map<int, std::vector<double>> IP_eff;
//    std::map<int, float> Total_eff_per_order;
//    std::vector<double> IP_total_eff;

//    std::map<int, tk::spline> tr_p;
    std::vector<tk::spline> tr_p;
//    std::map<int, tk::spline> tr_r;
    std::vector<tk::spline> tr_r;
//    std::map<int, tk::spline> tr_q;
    std::vector<tk::spline> tr_q;
//    std::map<int, tk::spline> tr_phi;
    std::vector<tk::spline> tr_phi;
//    std::map<int, tk::spline> tr_tx;
    std::vector<tk::spline> tr_tx;
//    std::map<int, tk::spline> tr_ty;
    std::vector<tk::spline> tr_ty;

//    std::map<int, std::vector<Eigen::Vector2f> > target_pos;

    PSF *psfs;
    Slit *slit;

    spectrograph_information spec_info;


};

#endif // MATRIXSIMULATOR_H
