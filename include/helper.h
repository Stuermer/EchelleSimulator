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

/*!
 * Saves vector to CSV File
 *
 * This function saves a std double vector to a CSV file e.g. for easy plotting with python
 * @param vec vector to be saved
 * @param filename path to output file
 */
void vectorToFile(std::vector<double> const& vec, std::string const& filename);

/*!
 * Saves an OpenCV matrix to a text file
 *
 * This function saves a cv::Mat matrix to a text file. The file is essentially a CSV file, but the matrix is 'numpy'
 * style, which means that is begins with [ and ends with ] characters.
 * @param image matrix to be saved
 * @param filename path to output file
 */
void MatToFile(cv::Mat& image, std::string const& filename);

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
std::vector<double> decompose_matrix(cv::Mat mat);

/*!
 * Composes 2x3 transformation matrix from shear, scale, rotation and translation parameters.
 *
 *
 * @param parameters [sx, sy, shear, \f$ \phi \f$, tx, ty]
 * @return 2x3 transformation matrix
 */
cv::Mat compose_matrix(std::vector<double> parameters);

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


double interpolate(const std::map<double,double> &data,
                    double x);

#endif