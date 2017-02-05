//
// Created by julian on 11.09.16.
//

#ifndef ECHELLESIMULATOR_NOISE_H
#define ECHELLESIMULATOR_NOISE_H
#include "opencv2/core.hpp"

class noise {
public:
    /**
     * Adds poisson noise (=photon noise) to the image.
     * @param image image where to add poisson noise
     */
    static void add_poisson(cv::Mat& image);

};


#endif //ECHELLESIMULATOR_NOISE_H
