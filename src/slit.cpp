#include "slit.h"
#include <opencv2/highgui.hpp>
#include <iostream>

slit::slit()
{
}

slit::slit(double w, double h, int slit_sampling){
    this->set_slit(w,h,slit_sampling);
}

void slit::set_slit(double w, double h, int slit_sampling){
    this->w = w;
    this->h = h;
    this->ratio = h/w;
    this->slit_sampling = slit_sampling;
    this->w_px = slit_sampling;
    this->h_px = slit_sampling * this->ratio;
    std::cout << (this->slit_image.type()) << std::endl ;

    #ifdef USE_GPU
    {
        cv::Mat ones =  cv::Mat::ones(round(this->h_px), this->w_px, CV_32FC1);
        this->slit_image = cv::gpu::GpuMat();
        this->slit_image.upload(ones);
    }
    #else
    {
        this->slit_image = cv::Mat::ones(round(this->h_px), this->w_px, CV_64F);
    }
    #endif
}

void slit::show(){
    cv::imshow("Slit Image", this->slit_image);
    cv::waitKey(0);

}
