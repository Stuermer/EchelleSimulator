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

struct raw_transformation{
    int order;
    double wavelength;
    cv::Mat transformation_matrix = cv::Mat(2,3,CV_64FC1);
    std::vector<double> decomposed_matrix;

};

struct spectrograph_information{
    double blaze;
    double gpmm;
};

typedef struct raw_transformation raw_transformation;

class MatrixSimulator
{
public:
    MatrixSimulator();
    void read_transformations(std::string path);
    cv::Mat get_transformation_matrix(int order, double wavelength);
    void calc_splines();
    void set_wavelength(int N);
    void set_wavelength(std::vector<double> wavelength);
    void calc_sim_matrices();
    #ifdef USE_GPU
    cv::gpu::GpuMat transform_slit(cv::gpu::GpuMat& slit_image, cv::Mat& transformation_matrix, double weight);
    int simulate_order(int order, cv::gpu::GpuMat& slit_image, cv::gpu::GpuMat& output_image, bool aberrations);
    void simulate_spectrum(cv::gpu::GpuMat& slit_image);
    #endif
    cv::Mat transform_slit(cv::Mat& slit_image, cv::Mat& transformation_matrix, double weight);
    int simulate_order(int order, cv::Mat& slit_image, cv::Mat& output_image, bool aberrations);
    void simulate_spectrum();

    void set_efficiencies(std::vector<Efficiency *> &efficiencies);
    void add_efficiency(Efficiency *eff);

    void add_source(Source * src);

    /**
     * Load spectrograph model from HDF file
     * @param path path to HDF file containing spectrograph model
     * @param fiber_number fiber to select
     */
    void load_spectrograph_model(std::string path, int fiber_number);

    
    void set_order_range(int min_order, int max_order);
    std::vector<int> orders;

    void set_ccd(CCD * ccd);
    void set_slit(Slit * slit);
    void set_psfs(PSF * psfs);

    int raw_n;


    void prepare_sources(std::vector<Source *> sources);
    void save_to_file(std::string filename, bool downsample=true, bool bleed=true, bool overwrite=false);
    void transformation_to_file(std::string filename);
    CCD * ccd;
    double get_blaze();
    double get_gpmm();
private:
    cv::Mat img;
    std::map<int, std::vector<raw_transformation> > raw_transformations;

    std::vector<Efficiency*> efficiencies;
    std::vector<Source*> sources;

    std::map<int, std::vector<double> > sim_wavelength;
    std::map<int, std::vector<cv::Mat> > sim_matrices;
    std::map<int, std::vector<double> > sim_efficiencies;
    std::map<int, std::vector<double> > sim_spectra;

    std::map<int, tk::spline > tr_p;
    std::map<int, tk::spline > tr_r;
    std::map<int, tk::spline > tr_q;
    std::map<int, tk::spline > tr_phi;
    std::map<int, tk::spline > tr_tx;
    std::map<int, tk::spline > tr_ty;

    PSF * psfs;
    Slit * slit;

    spectrograph_information spec_info;



};

#endif // MATRIXSIMULATOR_H
