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


inline double deg2rad(double deg)
{
    return deg * M_PI / 180.;
}

Efficiency::Efficiency()
{

}

Efficiency::~Efficiency()
{

}

std::vector< double > Efficiency::get_efficieny(int order, std::vector<double> wavelength)
{
    std::vector<double> res(wavelength.size(), 1.0);
    
    return res;

}


std::vector< double > GratingEfficiency::get_efficieny(int order, std::vector<double> wavelength)
{
  std::vector<double> res;  
  for(auto& w : wavelength)
    res.push_back(this->calc_eff(this->scalingfactor, order, deg2rad(this->alpha), deg2rad(this->blaze), w*1000., this->gpmm));
  return res;
}

double GratingEfficiency::calc_eff(double scalingfactor, int order, double alpha, double blaze, double wl, double n)
{
  double bb = asin(-sin(alpha)+(double) order*wl*n*1E-6);
  double x = (double) order * (cos(alpha) / cos(alpha - blaze)) * (cos(blaze) - sin(blaze) / tan((alpha + bb) / 2.));
  double sinc = sin(M_PI*x)/(M_PI*x);
  return scalingfactor * sinc*sinc;
                                       
}

GratingEfficiency::GratingEfficiency(double scalingfactor, double alpha, double blaze, double gpmm)
        : scalingfactor(scalingfactor), alpha(alpha), blaze(blaze), gpmm(gpmm) {}
