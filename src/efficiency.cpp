/*
 * <one line to give the library's name and an idea of what it does.>
 * Copyright (C) 2016  Julian Stürmer <julian.stuermer@online.de>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
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
                this->calc_eff(this->peak_efficiency, order, deg2rad(this->alpha), deg2rad(this->blaze), w, this->gpmm));
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

EtalonEfficiency::EtalonEfficiency(double d, double n, double theta, double R) : d(d / 1000.), n(n), theta(theta),
                                                                                 R(R) {
    this->cF = this->coefficient_of_finesse(R);
    this->integration_steps = 10;
}

EtalonEfficiency::~EtalonEfficiency() {

}

double EtalonEfficiency::coefficient_of_finesse(double R) {
    return 4. * R / ((1. - R) * (1. - R));
}

double EtalonEfficiency::T(double wl, double theta, double d, double n, double cF) {
    //delta = (2. * math.pi / wl) * 2. * n * math.cos(theta) * d
    //return 1. / (1. + cF * np.sin(0.5 * delta) ** 2)

    double delta = (2. * M_PI * n * cos(theta) * d) / wl;
    double sind = sin(delta);
    return 1. / (1. + cF * sind * sind);
}

double EtalonEfficiency::get_local_efficiency(double wavelength) {
    double res = this->T(wavelength / 1E6, theta, d, n, cF);
    return res;
}

std::vector<double> EtalonEfficiency::get_efficiency(int order, std::vector<double> &wavelength) {
    std::vector<double> spectrum;
    std::vector<double> diff;

    for (std::vector<int>::size_type i = 0; i != wavelength.size() - 1; i++) {
        diff.push_back(wavelength[i + 1] - wavelength[i]);
    }
    diff.push_back(diff.back());

    for (std::vector<int>::size_type i = 0; i != wavelength.size() - 1; i++) {
        double result = this->integral_s(wavelength[i] - diff[i] / 2., wavelength[i] + diff[i] / 2.,
                                         this->integration_steps);
        spectrum.push_back(result);
    }

    return spectrum;
}

std::vector<double> EtalonEfficiency::get_efficiency(int order, std::vector<double> &wavelength, int N) {

    std::vector<double> spectrum;
    std::vector<double> diff;

    for (std::vector<int>::size_type i = 0; i != N; i++) {
        diff.push_back(wavelength[i + 1] - wavelength[i]);
    }
    diff.push_back(diff.back());

    for (std::vector<int>::size_type i = 0; i != N; i++) {
        double result = this->integral_s(wavelength[i] - diff[i] / 2., wavelength[i] + diff[i] / 2.,
                                         this->integration_steps);
        spectrum.push_back(result);
    }

    return spectrum;
}

double EtalonEfficiency::integral_s(double a, double b, int n) {
    double step = (b - a) / n;  // width of each small rectangle
    float area = 0.0;  // signed area
    for (int i = 0; i < n; i++) {
        area += this->get_local_efficiency(a + (i + 0.5) * step) * step; // sum up each small rectangle
    }
    return area;
}

CSVEfficiency::CSVEfficiency(std::string path) {

    std::ifstream file(path.c_str());
    for (CSVIterator loop(file); loop != CSVIterator(); ++loop) {
        std::cout << (*loop)[0] << std::endl;
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