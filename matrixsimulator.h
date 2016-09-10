#ifndef MATRIXSIMULATOR_H
#define MATRIXSIMULATOR_H

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "spline.h"
#include "efficiency.h"

struct raw_transformation{
    int order;
    double wavelength;
    cv::Mat transformation_matrix = cv::Mat(2,3,CV_64FC1);
    std::vector<double> decomposed_matrix;

};

typedef struct raw_transformation raw_transformation;

class MatrixSimulator
{
public:
    MatrixSimulator();
    void read_transformations(std::string path);
    // QwtPlot * plot_transformations();
    cv::Mat get_transformation_matrix(int order, double wavelength);
    void calc_splines();
    void set_wavelength(int N);
    void calc_sim_matrices();
    cv::gpu::GpuMat transform_slit(cv::gpu::GpuMat& slit_image, cv::Mat& transformation_matrix, double weight);
    cv::Mat transform_slit(cv::Mat& slit_image, cv::Mat& transformation_matrix, double weight);
    int simulate_order(int order, cv::gpu::GpuMat& slit_image, cv::Mat& output_image);
    int simulate_order(int order, cv::Mat& slit_image, cv::Mat& output_image);
    cv::Mat simulate_spectrum(cv::gpu::GpuMat& slit_image);
    cv::Mat simulate_spectrum(cv::Mat& slit_image);
    
    void prepare_efficienies(std::vector<GratingEfficiency>& efficiencies);
    
    

//private:
    std::map<int, std::vector<raw_transformation> > raw_transformations;
    std::vector<int> orders;

    std::map<int, std::vector<double> > sim_wavelength;
    std::map<int, std::vector<cv::Mat> > sim_matrices;
    std::map<int, std::vector<double> > sim_efficiencies;


    std::map<int, tk::spline > tr_p;
    std::map<int, tk::spline > tr_r;
    std::map<int, tk::spline > tr_q;
    std::map<int, tk::spline > tr_phi;
    std::map<int, tk::spline > tr_tx;
    std::map<int, tk::spline > tr_ty;

    int raw_n;

};

#endif // MATRIXSIMULATOR_H
