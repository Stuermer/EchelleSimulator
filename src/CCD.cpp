//
// Created by julian on 21.09.16.
//

#include <opencv2/imgproc.hpp>
#include "CCD.h"
#include "hdf5opencv.h"
#include <CCfits/CCfits>

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

void CCD::save_to_hdf(std::string filename, bool downsample, bool bleed, bool overwrite) {
    cv::Mat res = this->get_image(downsample, bleed);
    hdf5opencv::hdf5save(filename.c_str(), "image", res, overwrite);
    // MatToFile(res, filename);
}

void CCD::save_to_fits(std::string filename, bool downsample, bool bleed, bool overwrite) {
    cv::Mat res = this->get_image(downsample, bleed);
    long naxis    =   2;
    long naxes[2] = { res.cols, res.rows };
    std::auto_ptr<CCfits::FITS> pFits(0);
    if (overwrite){
        remove(filename.c_str());
    }
    try
    {
        // overwrite existing file if the file already exists.

        const std::string fileName(filename.c_str());

        // Create a new FITS object, specifying the data type and axes for the primary
        // image. Simultaneously create the corresponding file.

        // this image is unsigned short data, demonstrating the cfitsio extension
        // to the FITS standard.

        pFits.reset( new CCfits::FITS(filename, FLOAT_IMG , naxis , naxes ) );
    }
    catch (CCfits::FITS::CantCreate)
    {
        // ... or not, as the case may be.
        std::cout<<"Can't create FITS file."<< std::endl;
    }
    long nelements(1);
    nelements = std::accumulate(&naxes[0],&naxes[naxis],1,std::multiplies<long>());
    long  fpixel(1);
    std::valarray<double> data_array(nelements);
    for (int i=0; i<res.rows; ++i){
        for (int j=0; j<res.cols; ++j){
            data_array[j*res.cols+i] = res.at<double>(j,i);
        }

    }

    pFits->pHDU().write(fpixel, nelements, data_array);
    pFits->flush();

}
CCD::~CCD() {

}

cv::Mat CCD::get_image(bool downsample, bool bleed) {
    cv::Mat result;
    if (downsample){
        #ifdef  USE_GPU
                {
                this->data.download(result)
            };
        #else
                {
                    result = cv::Mat(Ny, Nx, this->data.type());
                }
        #endif
        cv::resize(this->data, result, result.size(), cv::INTER_NEAREST);
    }
    else {
        result = this->data;
    }

    if (bleed){
        this->do_bleed(result, 650000.);
    }

    return result;
}

void do_bleed_updown(cv::Mat &input, double & limit, int i, int j){
    //find limit
    int k = i;
    while ((k>1) && (k<input.rows-1) && (input.at<double>(k, j)))
    {
        k++;
    }

    cv::blur(input.colRange(i-1,k+1), input.colRange(i-1,k+1), cv::Size(1,3),cv::Point(0.,1.5) ,cv::BORDER_REPLICATE);
    if(i>0)
    {
        if (input.at<double>(i-1,j) > limit)
            do_bleed_updown(input, limit, i-1, j);
    }
//    double diff = input.at<double>(i, j) - limit;
//    if ((diff > 0) && (i + 2 < input.rows) && (i - 2 > 0)) {
//        std::cout<< i << "\t" << j << std::endl;
//        input.at<double>(i + 1, j) += diff/2.;
//        input.at<double>(i - 1, j) += diff/2.;
////        input.at<double>(i - 1, j) += diff / 2.;
//        input.at<double>(i, j) -= diff;
//        do_bleed_updown(input, limit, i - 1, j);
//        do_bleed_updown(input, limit, i + 1, j);
//    }
}

//void do_bleed_down(cv::Mat &input, double limit, int i, int j){
//    double diff = input.at<double>(i,j) - limit;
//    if ((diff>0) && (i+1 < input.rows) && (i-1 >0))  {
//        input.at<double>(i - 1, j) += diff / 2.;
//        do_bleed_down(input, limit, i - 1, j);
//    }
//
//}

void CCD::do_bleed(cv::Mat &input, double limit) {

    for(int i=0; i<input.cols; ++i) {
        for (int j = 0; j < input.rows; ++j) {
            if(input.at<double>(j,i) > limit){
                do_bleed_updown(input, limit, j, i);
            }
        }
    }
}

