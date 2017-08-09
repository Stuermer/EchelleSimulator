//
// Created by zrobertson on 8/8/17.
//

#include "telescope.h"
#include <math.h>

Telescope::Telescope() {
    this->d = 2/sqrt(M_PI);
    this->f = 1.;
}

Telescope::~Telescope() {

}

Telescope::Telescope(double d, double f): d(d),f(f) {

}

double Telescope::get_area() {
    return this->d*this->d / 4. * M_PI;
}

double Telescope::get_diameter() {
    return this->d;
}

double Telescope::get_focal_ratio() {
    return this->f / this->d;
}