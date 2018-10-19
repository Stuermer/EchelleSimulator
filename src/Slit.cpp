#include "Slit.h"
#include <iostream>

Slit::Slit()
= default;

Slit::Slit(double w, double h, int slit_sampling) {
    this->set_slit(w, h, slit_sampling);
}

void Slit::set_slit(double w, double h, int slit_sampling) {
    this->w = w;
    this->h = h;
    this->ratio = h / w;
    this->slit_sampling = slit_sampling;
    this->w_px = slit_sampling;
    this->h_px = static_cast<int>(slit_sampling * this->ratio);
    this->slit_image = Matrix(static_cast<size_t>(this->h_px), static_cast<size_t>(this->w_px));
}

