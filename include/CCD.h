#ifndef ECHELLESIMULATOR_CCD_H
#define ECHELLESIMULATOR_CCD_H

#include "hdf5.h"
#include <string>
#include <vector>

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

    /**
     * Saves simulated data to HDF file
     * @param filename file path
     * @param bleed bleed overexposed pixel TODO: not implemented yet
     * @param overwrite set true to overwrite existing file
     */
    void save_to_hdf(std::string filename, bool bleed = true, bool overwrite = false);

    /**
     * Saves simulated data to FITS file
     * @param filename file path
     * @param overwrite set true to overwrite existing file
     */
    void save_to_fits(std::string filename, bool overwrite);

    /**
     * Returns pixel of the CCD detector size in microns.
     * @return Pixel size [microns]
     */
    double *get_pixelsize();

    //overload + operator
//    static void do_bleed(cv::Mat &input, double limit);

    std::vector<int> data;
    bool use_gpu = false;
    int Nx, Ny;
    double pixel_size;


};


#endif //ECHELLESIMULATOR_CCD_H
