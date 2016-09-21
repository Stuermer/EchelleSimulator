//
// Created by julian on 21.09.16.
//

#ifndef ECHELLESIMULATOR_CCD_H
#define ECHELLESIMULATOR_CCD_H


#include <opencv2/core/types_c.h>
#include "hdf5.h"
#include <string>
#include <ml.h>

/*!
 * \class CCD glass
 * \brief class representing a CCD detector
 */
class CCD {
public:
    /*!
     * Constructor
     * @param Nx number of pixels in X direction
     * @param Ny number of pixels in Y direction
     * @param data_type data type, should be same as slit.image. Possible values see opencv datatypes (e.g. CV_32F, CV_64F,...)
     * @return CCD
     */
    CCD(int Nx, int Ny, int oversampling, int data_type);
    ~CCD();
    void save_to_file(std::string filename, bool downsample=true, bool overwrite=false);

    cv::Mat get_image(bool downsample=true);
    //overload + operator
    CCD operator+(const CCD & ccd){
    this->data += ccd.data;
    }
    cv::Mat data;


private:
    int Nx, Ny, oversampling;



};


#endif //ECHELLESIMULATOR_CCD_H
