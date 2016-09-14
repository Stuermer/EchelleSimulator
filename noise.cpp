//
// Created by julian on 11.09.16.
//

#include "noise.h"
#include <random>

void noise::add_poisson(cv::Mat &image) {

    std::default_random_engine generator;
    for(int i=0; i<image.rows; ++i){
        for(int j=0; j<image.cols; ++j)
        {
            if (image.at<float>(j,i) > 1E-10){
                std::poisson_distribution<int> distribution(image.at<float>(j,i));
                int value = distribution(generator);
                image.at<float>(j, i) = value;

            }
        }
    }




}
