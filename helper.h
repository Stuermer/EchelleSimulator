#ifndef HELPER_H
#define HELPER_H
#define _USE_MATH_DEFINES

#include <vector>
#include <cmath>
#include "opencv2/core.hpp"

std::vector<double> decompose_matrix(cv::Mat mat){
    std::vector<double> result;

    double a = mat.at<double>(0,0);
    double b = mat.at<double>(0,1);
    double d = mat.at<double>(1,0);
    double e = mat.at<double>(1,1);
    double tx = mat.at<double>(0,2);
    double ty = mat.at<double>(1,2);
    // p=sqrt(a^2+b^2)
    double p = sqrt(a*a+b*b);
    result.push_back(p);
    // r = ae-bd / p
    result.push_back((a*e-b*d)/p);
    result.push_back((a*d+b*e)/(a*e-b*d));
    double phi = atan2(b,a);
    if (phi>0)
      phi -= 2.*M_PI;
    result.push_back(phi);
    result.push_back(tx);
    result.push_back(ty);
    // return <p,r,q,phi, tx, ty>
    return result;
}

cv::Mat compose_matrix(std::vector<double> parameters){
    cv::Mat a = cv::Mat::zeros(2,2,CV_64FC1);
    cv::Mat tmp1 = cv::Mat::zeros(2,2,CV_64FC1);
    cv::Mat tmp2 = cv::Mat::zeros(2,2,CV_64FC1);
    cv::Mat tmp3 = cv::Mat::zeros(2,2,CV_64FC1);

    tmp1.at<double>(0,0) = parameters[0];
    tmp1.at<double>(1,1) = parameters[1];

    tmp2.at<double>(0,0) = 1.;
    tmp2.at<double>(1,1) = 1.;
    tmp2.at<double>(1,0) = parameters[2];

    double cosp = cos(parameters[3]);
    double sinp = sin(parameters[3]);

    tmp3.at<double>(0,0) = cosp;
    tmp3.at<double>(0,1) = sinp;
    tmp3.at<double>(1,0) = -sinp;
    tmp3.at<double>(1,1) = cosp;

    a = tmp1 * tmp2 * tmp3;

    cv::Mat result = cv::Mat(2,3,CV_64FC1);

    result.at<double>(0,0) = a.at<double>(0,0);
    result.at<double>(0,1) = a.at<double>(0,1);
    result.at<double>(1,0) = a.at<double>(1,0);
    result.at<double>(1,1) = a.at<double>(1,1);
    result.at<double>(0,2) = parameters[4];
    result.at<double>(1,2) = parameters[5];

    return result;
}

#endif // HELPER_H
