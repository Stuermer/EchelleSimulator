#ifndef HELPER_H
#define HELPER_H
/// load pi
#define _USE_MATH_DEFINES

#include <vector>
#include <cmath>
#include <vector>
#include <iterator>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cstdlib>
#include "csv_reader.h"
#include "H5Cpp.h"
#include <map>
#include <array>
#include <algorithm>
#include <numeric>
#include <random>

/**
 * Saves vector to CSV File
 *
 * This function saves a std double vector to a CSV file e.g. for easy plotting with python
 * @param vec vector to be saved
 * @param filename path to output file
 */
void vector_to_file(std::vector<double> const &vec, std::string const &filename);

/**
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
std::array<double, 6> decompose_matrix(std::array<double, 6> mat);

/**
 * Composes 2x3 transformation matrix from shear, scale, rotation and translation parameters.
 *
 *
 * @param parameters [sx, sy, shear, \f$ \phi \f$, tx, ty]
 * @return 2x3 transformation matrix
 */
std::array<double, 6> compose_matrix(std::vector<double> parameters);

/**
 * Calculates sorted index array of a given vector.
 *
 * Equivalent to numpy.argsort
 *
 * @param v vector to be sorted
 * @return array of indices that would sort the vector v
 */
std::vector<std::size_t> compute_sort_order(const std::vector<double> &v);

/**
 * Wraps an angle around [-pi, pi]
 * @param r in radian
 * @return wrapped angle
 */
double wrap_rads(double r);

/**
 * Linear interpolation of data at position x
 * @param data data to interpolate
 * @param x x
 * @return interpolated value
 */
double interpolate(const std::map<double, double> &data, double x);

/*
 * Helper for HDF file to get group names
 */
herr_t file_info(hid_t loc_id, const char *name, const H5L_info_t *linfo, void *opdata);


/**
 * 
 * @tparam Out 
 * @param s 
 * @param delim 
 * @param result 
 */
template<typename Out>
void split(const std::string &s, char delim, Out result);

/**
 * Split string at delimiter and fill vector
 * @param s string to split
 * @param delim delimiter
 * @return vector of strings
 */
std::vector<std::string> split_to_vector(const std::string &s, char delim);

/**
 * Downloads phoenix spectrum for given parameters.
 *
 * @param t_eff effective temperature [Kelvin]
 * @param log_g surface gravity
 * @param z overall metallicity
 * @param alpha alpha element abundance
 * @param path where to save the spectrum
 * @return 0 if download succeeded, CURLeCode otherwise
 */
int download_phoenix(int t_eff, double log_g, double z, double alpha, const std::string path);

/**
 * Downloads phoenix wavelength grid file
 * @param path where to save wavelength grid file
 * @return 0 if succeeded, CURLeCode otherwise
 */
int download_wave_grid(std::string path);

/**
 * Checks if file exists
 * @param path file path to check
 * @return true/false if file exists/doesn't exist
 */
bool check_for_file(const std::string &path);

/*
 * Helper for curl to write data to hard drive
 */
size_t write_data(void *ptr, size_t size, size_t nmemb, FILE *stream);

#endif
