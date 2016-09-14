//
// Created by julian on 11.09.16.
//

#ifndef ECHELLESIMULATOR_NOISE_H
#define ECHELLESIMULATOR_NOISE_H
#include "opencv2/core.hpp"

class noise {
public:
    static void add_poisson(cv::Mat& image);

};


#endif //ECHELLESIMULATOR_NOISE_H
