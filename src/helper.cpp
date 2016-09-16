#include "helper.h"

#include <vector>
#include <cmath>
#include "opencv2/core.hpp"

#include <vector>
#include <iterator>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <opencv2/imgproc.hpp>
#include <highgui.h>


void vectorToFile(std::vector<double> const& vec, std::string const& filename) {
  std::ofstream file(filename);
  auto first = true;
  for (float f : vec) 
  { 
      if (!first) { file << ","; } 
      first = false; 
      file << f; 
  }
  file << std::endl;
  file.close();  
}

void MatToFile(cv::Mat& image, std::string const& filename){
  std::ofstream file(filename);
  
  file << cv::format(image, "numpy") << std::endl;
  file.close();
}

std::vector<double> decompose_matrix(cv::Mat mat){
    std::vector<double> result;
    /*
     * Matrix looks like:
     * a b tx
     * d e ty
     */

    double a = mat.at<double>(0,0);
    double b = mat.at<double>(0,1);
    double d = mat.at<double>(1,0);
    double e = mat.at<double>(1,1);
    double tx = mat.at<double>(0,2);
    double ty = mat.at<double>(1,2);
    
    double sx = sqrt(a*a+d*d);
    double sy = sqrt(b*b+e*e);
    
    double phi = atan2(b,a);
    if (phi>0)
      phi -= 2.*M_PI;	
    
    double shear = atan2(-d,e) - phi;
    result.push_back(sx);
    result.push_back(sy);
    result.push_back(shear);
    result.push_back(phi);
    result.push_back(tx);
    result.push_back(ty);
    // return <sx, sy, shear, rot, tx ,ty>
    // return <p,r,q,phi, tx, ty>
    return result;
}

cv::Mat compose_matrix(std::vector<double> parameters){
  double sx = parameters[0];
  double sy = parameters[1];
  double shear = parameters[2];
  double rot = parameters[3];

  cv::Mat result = cv::Mat(2,3,CV_64FC1);
    
  result.at<double>(0,0) = sx*cos(rot);
  result.at<double>(0,1) = -sy*sin(rot+shear);
  result.at<double>(1,0) = sx*sin(rot);
  result.at<double>(1,1) = sy*cos(rot+shear);
  result.at<double>(0,2) = parameters[4];
  result.at<double>(1,2) = parameters[5];

  return result;
}

std::vector<std::size_t> compute_order(const std::vector<double> &v) {
    std::vector<std::size_t> indices(v.size());
    std::iota(indices.begin(), indices.end(), 0u);
    std::sort(indices.begin(), indices.end(), [&](int lhs, int rhs) {
        return v[lhs] < v[rhs];
    });
    std::vector<std::size_t> res(v.size());
    for (std::size_t i = 0; i != indices.size(); ++i) {
        res[indices[i]] = i;
    }
    return res;
}

void show_cv_matrix(cv::Mat img, std::string windowname="image") {
    double minVal, maxVal;

    cv::Mat img_show = img.clone();
    cv::minMaxLoc(img_show, &minVal, &maxVal); //find minimum and maximum intensities

    cv::cvtColor(img_show, img_show, CV_GRAY2RGB);
    cv::namedWindow(windowname,CV_WINDOW_NORMAL);
    cv::imshow(windowname, img_show);
    cv::waitKey(0);

}

void print_cv_matrix_info(cv::Mat img, std::string imagename="Image") {
    std::cout<< "------------" << imagename << "----------------:" << std::endl;
    std::cout << "Dimensions: " << img.rows << " x " << img.cols << std::endl;
    std::cout << "Type: " << img.type() << std::endl;
    double minVal, maxVal;
    cv::minMaxLoc(img, &minVal, &maxVal); //find minimum and maximum intensities
    std::cout << "Min/Max: " << minVal << " / " << maxVal << std::endl;
    std::cout<< "------------------------------------------:"<<std::endl << std::endl;


}

