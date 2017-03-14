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

#ifndef EFFICIENCY_H
#define EFFICIENCY_H
#include <string>
#include <vector> 


class Efficiency
{
public:
    Efficiency();
    virtual ~Efficiency();
    virtual std::vector<double> get_efficieny(int order, std::vector<double> wavelength);
private:
  
};
class ConstantEfficiency: public Efficiency
{
public:
    ConstantEfficiency(double efficiency);
    std::vector<double> get_efficienct(int order, std::vector<double> wavelength);
private:
    double eff;
};

class GratingEfficiency : public Efficiency
{
public:
    GratingEfficiency(double scalingfactor, double alpha, double blaze, double gpmm);

    std::vector<double> get_efficieny(int order, std::vector<double> wavelength);

private:
  double scalingfactor=0.8;
  double alpha=76.;
  double blaze=76.;
  double gpmm=31.6;
  double calc_eff(double scalingfactor, int order, double alpha, double blaze, double wl, double n);
      
};

#endif // EFFICIENCY_H
