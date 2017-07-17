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
#include "csv_reader.h"
#include "H5Cpp.h"
#include "Eigen/Dense"
#include <map>

typedef Eigen::Matrix<float, 2, 3> Matrix23f;

/*!
 * Saves vector to CSV File
 *
 * This function saves a std double vector to a CSV file e.g. for easy plotting with python
 * @param vec vector to be saved
 * @param filename path to output file
 */
void vectorToFile(std::vector<double> const &vec, std::string const &filename);

/*!
 * Decomposes a 2x3 affine transformation matrix into its underlying geometric components.
 *
 * A transformation matrix \f$ \begin{pmatrix} a & b & c \\ d & e & f \end{pmatrix} \f$ can be decomposed in
 * rotation \f$ \phi \f$, scale in X- and Y direction \f$ sx, sy \f$, shearing \f$ s \f$ and translation in X- and Y
 * \f$ tx, ty \f$.
 *
 * The decomposition is not unique, there exist other decompositions. However, it seems that this particular decomposition
 * is rather stable and for each component one gets a rather smooth function across an order.
 *
 * \warning When rotation \f$ \phi \f$ is close to 0, it sometimes alters between \f$ \pm 2\pi \f$. This is a problem when
 * trying to make a smooth spline across an order. Temporary fix is to add \f$ 2 \pi \f$ when values are near \f$ -\pi \f$.
 *
 * @param mat 2x3 transformation matrix
 * @return [sx, sy, shear, \f$ \phi \f$, tx, ty]
 */
std::vector<double> decompose_matrix(Matrix23f mat);

/*!
 * Composes 2x3 transformation matrix from shear, scale, rotation and translation parameters.
 *
 *
 * @param parameters [sx, sy, shear, \f$ \phi \f$, tx, ty]
 * @return 2x3 transformation matrix
 */
Matrix23f compose_matrix(std::vector<double> parameters);

/*!
 * Calculates sorted index array of a given vector.
 *
 * Equivalent to numpy.argsort
 *
 * @param v vector to be sorted
 * @return array of indices that would sort the vector v
 */
std::vector<std::size_t> compute_sort_order(const std::vector<double> &v);

/*!
 * Plots an OpenCV matrix of type CV_32.
 *
 * For plotting a double/float matrix, the matrix needs to be rescaled before plotted on screen.
 * @param img Matrix/Image to be plotted
 * @param windowname name of the window containing the plot
 */
void show_cv_matrix(cv::Mat img, std::string windowname);

/*!
 * Prints out basic information about a Matrix/Image.
 *
 * Prints type, dimensions, min- and max value.
 *
 * @param img Matrix/Image to be evaluated.
 * @param imagename Name of the matrix
 */
void print_cv_matrix_info(cv::Mat img, std::string imagename);

/**
 * Wraps an angle around [-pi, pi]
 * @param r in radian
 * @return wrapped angle
 */
double wrap_rads(double r);

/**
 * Creates fits file
 * @param filename
 */
void create_fits_file(std::string filename);

double interpolate(const std::map<double, double> &data, double x);

herr_t file_info(hid_t loc_id, const char *name, const H5L_info_t *linfo, void *opdata);


std::vector<float> random_from_2_distributions(std::vector<float> wl, std::vector<float> density1, std::vector<float> density2, int N_samples);

//int save_to_fits(const std::string filename, cv::Mat img);
#endif
