#ifndef HELPER_H
#define HELPER_H
#define _USE_MATH_DEFINES

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


void vectorToFile(std::vector<double> const& vec, std::string const& filename);
void MatToFile(cv::Mat& image, std::string const& filename);

std::vector<double> decompose_matrix(cv::Mat mat);

cv::Mat compose_matrix(std::vector<double> parameters);

std::vector<std::size_t> compute_order(const std::vector<double>& v);

void show_cv_matrix(cv::Mat img, std::string windowname);

void print_cv_matrix_info(cv::Mat img, std::string imagename);
#endif // HELPER_H
