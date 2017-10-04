//
// Created by julian on 21.09.16.
//

#ifndef ECHELLESIMULATOR_CCD_H
#define ECHELLESIMULATOR_CCD_H


#include "opencv2/core/types_c.h"
#include "hdf5.h"
#include <string>
#include <ml.h>

#ifdef USE_GPU
#include "opencv2/gpu/gpu.hpp"
#endif

/*!
 * \class CCD
 * \brief class representing a CCD detector
 */
class CCD {
public:
    /*!
     * Constructor
     * @param Nx number of pixels in X direction
     * @param Ny number of pixels in Y direction
     * @param pixelsize size of the pixels [micron]
     * @return CCD
     */
    CCD(int Nx, int Ny, double pixelsize);

    ~CCD();

    void save_to_hdf(std::string filename, bool downsample = true, bool bleed = true, bool overwrite = false);

    void save_to_fits(std::string filename, bool downsample = true, bool bleed = true, bool overwrite = false);

    cv::Mat get_image(bool downsample = true, bool bleed = true);

    double * get_pixelsize();

    //overload + operator
    static void do_bleed(cv::Mat &input, double limit);

#ifdef USE_GPU
    CCD operator+(const CCD & ccd){
        cv::gpu::add(this->data, ccd.data, this->data);
    }
#else

    CCD operator+(const CCD &ccd) {
        this->data += ccd.data;
    }

#endif

#ifdef USE_GPU
    cv::gpu::GpuMat data;
    bool use_gpu = true;
#else
    cv::Mat data;
    bool use_gpu = false;
#endif


private:
    int Nx, Ny;
    double pixelsize;


};


#endif //ECHELLESIMULATOR_CCD_H
