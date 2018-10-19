#define _USE_MATH_DEFINES

#include "efficiency.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <csv_reader.h>
#include <helper.h>


inline double deg2rad(double deg) {
    return deg * M_PI / 180.;
}

Efficiency::Efficiency() {

}

Efficiency::~Efficiency() {

}

std::vector<double> Efficiency::get_efficiency(int order, std::vector<double> &wavelength) {
    std::vector<double> res(wavelength.size(), 1.0);
    return res;
}

std::vector<double> Efficiency::get_efficiency(int order, std::vector<double> &wavelength, int N) {
    std::vector<double> res(N, 1.0);
    return res;
}


std::vector<double> GratingEfficiency::get_efficiency(int order, std::vector<double> &wavelength) {
    std::vector<double> res;
    for (auto &w : wavelength)
        res.push_back(
                this->calc_eff(this->peak_efficiency, order, deg2rad(this->alpha), deg2rad(this->blaze), w,
                               this->gpmm));
    return res;
}

std::vector<double> GratingEfficiency::get_efficiency(int order, std::vector<double> &wavelength, int N) {
    std::vector<double> res;
    for (std::vector<int>::size_type i = 0; i != N; i++)
        res.push_back(
                this->calc_eff(this->peak_efficiency, order, deg2rad(this->alpha), deg2rad(this->blaze), wavelength[i],
                               this->gpmm));
    return res;
}

double
GratingEfficiency::calc_eff(double scalingfactor, int order, double alpha, double blaze, double wl, double gpmm) {
    double bb = asin(-sin(alpha) + (double) order * wl * 1E-6 / (1. / gpmm / 1000.));
    double x = (double) order * (cos(alpha) / cos(alpha - blaze)) * (cos(blaze) - sin(blaze) / tan((alpha + bb) / 2.));
    double sinc = sin(M_PI * x) / (M_PI * x);
    return scalingfactor * sinc * sinc;

}

GratingEfficiency::GratingEfficiency(double peak_efficiency, double alpha, double blaze, double gpmm)
        : peak_efficiency(peak_efficiency), alpha(alpha), blaze(blaze), gpmm(gpmm) {}

ConstantEfficiency::ConstantEfficiency(double efficiency) : eff(efficiency) {

}

std::vector<double> ConstantEfficiency::get_efficiency(int order, std::vector<double> &wavelength) {
    std::vector<double> res;
    for (auto &w : wavelength)
        res.push_back(this->eff);
    return res;
}

std::vector<double> ConstantEfficiency::get_efficiency(int order, std::vector<double> &wavelength, int N) {
    std::vector<double> res;
    for (std::vector<int>::size_type i = 0; i != N; i++)
        res.push_back(this->eff);
    return res;
}


CSVEfficiency::CSVEfficiency(std::string path) {

    std::ifstream file(path.c_str());
    for (CSVIterator loop(file); loop != CSVIterator(); ++loop) {
        wl.push_back(stod((*loop)[0]));
        ef.push_back(stod((*loop)[1]));
        this->data.insert(std::pair<double, double>(stod((*loop)[0]), stod((*loop)[1])));
    }

}

std::vector<double> CSVEfficiency::get_efficiency(int order, std::vector<double> &wavelength) {
    std::vector<double> res;
    for (auto &w : wavelength)
        res.push_back(interpolate(this->data, w));
    return res;
};

std::vector<double> CSVEfficiency::get_efficiency(int order, std::vector<double> &wavelength, int N) {
    std::vector<double> res;
    for (std::vector<int>::size_type i = 0; i != N; i++) {
        res.push_back(interpolate(this->data, wavelength[i]));
    }
    return res;
};