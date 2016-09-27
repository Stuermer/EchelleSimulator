//
// Created by julian on 21.09.16.
//

#include <opencv2/imgproc.hpp>
#include "CCD.h"
#include "hdf5opencv.h"

CCD::CCD(int Nx, int Ny, int oversampling, int data_type) :Nx(Nx), Ny(Ny), oversampling(oversampling) {

#ifdef USE_GPU
    {
        cv::Mat ones =  cv::Mat::zeros(Ny*oversampling, Nx*oversampling, CV_32FC1);
        this->data = cv::gpu::GpuMat();
        this->data.upload(ones);
    }
#else
    {
        this->data = cv::Mat::zeros(Ny*oversampling, Nx*oversampling, CV_64F);
    }
#endif

    // this->data = cv::Mat::zeros(Ny*oversampling, Nx*oversampling, data_type);

}

void CCD::save_to_file(std::string filename, bool downsample, bool overwrite) {
    cv::Mat res = this->get_image(downsample);
    hdf5opencv::hdf5save(filename.c_str(), "image", res, overwrite);
    // MatToFile(res, filename);
}

CCD::~CCD() {

}

cv::Mat CCD::get_image(bool downsample) {
    if (downsample){
        cv::Mat result = cv::Mat(Ny, Nx, this->data.type());
        cv::resize(this->data, result, result.size(), cv::INTER_NEAREST);
        return result;
    }
    else {
#ifdef USE_GPU
    {
        cv::Mat result;
        this->data.download(result);
        return result;
    }
#else
        {
            return this->data;
        }
#endif
}
}