//
// Created by julian on 13.09.16.
//

#ifndef ECHELLESIMULATOR_PSF_H
#define ECHELLESIMULATOR_PSF_H

#include <string>
#include <map>
#include "opencv2/core.hpp"

struct PSFdata
{
    double wavelength;
    cv::Mat psf;

    PSFdata(double w, cv::Mat p) : wavelength(w), psf(p.clone()) {};

    bool operator < (const PSFdata& str) const
    {
        return (wavelength < str.wavelength);
    }
};

/*!
 * \class PSF
 * \brief Class handles point spread functions
 *
 * This class handles point spread functions (PSFs). It's basic functionality is to deliver a PSF as a 2d matrix for a
 * given wavelength and order of an echelle spectrograph.
 *
 * Typically the PSF of an echelle spectrograph will vary across the CCD depending on the wavelength, the echelle order,
 * and the illumination of the optics. The PSF might not be stable from target to target, as illumination might vary due
 * to different coupling conditions and imperfect scrambling of the fibers.
 *
 * To implement own PSF classes, inherit from PSF and overwrite get_PSF function
 */
class PSF {
public:
    PSF();
    virtual ~PSF();
    virtual cv::Mat get_PSF(int order, double wavelength) = 0;
};

/**
 * \class PSF_ZEMAX
 * \brief class representing a PSF as computed by ZEMAX
 *
 * Zemax returns a 2D array with a calculated Huygens PSF. It is usually saved in the spectrograph model.
 */
class PSF_ZEMAX : public PSF{
public:
    /**
     * Creator
     * @param filename filename of the spectrograph model
     * @param fiber_number fiber number for the PSF model
     */
    PSF_ZEMAX(std::string filename, int fiber_number);
    cv::Mat get_PSF(int order, double wavelength);
private:
    /**
     * Interpolates between neighboring PSFs as a simple linear sum
     * @param psf1 first PSF
     * @param psf2 second PSF
     * @param w1 first PSF wavelength
     * @param w2 second PSF wavelength
     * @param w wavelength to use for calculating PSF
     * @return
     */
    cv::Mat interpolate_PSF(cv::Mat psf1, cv::Mat psf2, double w1, double w2, double w);
    std::map< int, std::vector<PSFdata> > psfs;

};

/**
 * \class PSF_gaussian
 * \brief Class representing a gaussian PSF
 *
 * This class handels gaussian like PSF functions.
 */
class PSF_gaussian : public PSF{
public:
    /**
     * Creator
     * @param sigma Sigma of the gaussian in (oversampled) pixel
     * @param aperture size of the gaussian kernel (in number of oversampled pixels)
     */
    PSF_gaussian(double sigma, double aperture=3.);

    /**
     * Returns gaussian PSF, independently of order and wavelength
     * @param order echelle order
     * @param wavelength  wavelength in microns
     * @return
     */
    cv::Mat get_PSF(int order, double wavelength);

private:
    double sigma;
    int ksize;

};

#endif //ECHELLESIMULATOR_PSF_H
