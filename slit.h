#ifndef SLIT_H
#define SLIT_H
#include <opencv2/core.hpp>
#include "opencv2/gpu/gpu.hpp"

class slit
{
public:
    slit();
    slit(double w, double h, int slit_sampling);
    void set_slit(double w, double h, int slit_sampling);
    void show();

//private:
    double w;
    double h;
    int w_px;
    int h_px;
    double ratio;

    int slit_sampling;

    #ifdef USE_GPU
        cv::gpu::GpuMat slit_image;
        bool use_gpu = true;
    #else
        cv::Mat slit_image;
        bool use_gpu = false;
    #endif
};

#endif // SLIT_H
