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

#ifndef SOURCE_H
#define SOURCE_H

#include <vector>
#include <map>

class Source
{
public:
    Source();
    virtual ~Source();
    virtual std::vector<double> get_spectral_density(std::vector<double> wavelength);
    std::vector<double> get_spectrum(std::vector<double> wavelength);
    virtual double get_spectral_density(double wavelength);


private:
    double integral_s(double a, double b, int n);

};

class Constant : public Source{
public:
    Constant();
    Constant(double value);
    std::vector<double> get_spectral_density(std::vector<double> wavelength);
    double get_spectral_density(double wavelength);
private:
    double value;
};

class IdealEtalon : public Source{
public:
    IdealEtalon(double d, double n, double theta, double R);
    double coefficient_of_finesse();
    double FSR();
    double F();
    static double T(double wl, double theta, double d, double n, double cF);
    std::vector<double> get_spectral_density(std::vector<double> wavelength);
    double get_spectral_density(double wavelength);
private:
    double d;
    double n;
    double theta;
    double R;
    double cF;
};
#endif // SOURCE_H
