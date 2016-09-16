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


class PSF {
public:
    PSF(std::string filename);
    cv::Mat get_PSF(int order, double wavelength);

    cv::Mat interpolate_PSF(cv::Mat psf1, cv::Mat psf2, double w1, double w2, double w);
private:
    std::map< int, std::vector<PSFdata> > psfs;
};


#endif //ECHELLESIMULATOR_PSF_H
