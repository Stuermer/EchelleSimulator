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

#include <cmath>
#include <algorithm>
#include "source.h"
#include <CCfits/CCfits>
#include <helper.h>

Source::Source()
{
    //default value for number of integration steps.
    this->integration_steps = 10;
    this->shift = 1.;
}

Source::~Source()
{

}

std::vector<double> Source::get_spectral_density(std::vector<double> wavelength) {
    return std::vector<double>();
}

double Source::get_spectral_density(double wavelength) {
    return 1.;
}
std::vector<double> Source::get_spectrum(std::vector<double> wavelength) {
    std::vector<double> spectrum;
    std::vector<double> diff;

    for(std::vector<int>::size_type i = 0; i != wavelength.size()-1; i++)
    {
        diff.push_back(this->shift*(wavelength[i+1] - wavelength[i]));
    }
    diff.push_back(diff.back());

    for(std::vector<int>::size_type i = 0; i != wavelength.size()-1; i++)
    {
        double result = this->integral_s( this->shift*wavelength[i] - diff[i] / 2., this->shift*wavelength[i] + diff[i] / 2., this->integration_steps);
        spectrum.push_back(result);
    }

    return spectrum;
}

void Source::set_integration_steps(int n) {
    this->integration_steps = n;
}

double Source::integral_s(double a, double b, int n) {
        double step = (b - a) / n;  // width of each small rectangle
        double area = 0.0;  // signed area
        for (int i = 0; i < n; i ++) {
            area += 10000000.* this->get_spectral_density(a + (i + 0.5) * step) * step; // sum up each small rectangle
        }
        return area;
}

void Source::set_doppler_shift(double shift) {
    this->shift = 1. + shift / 300000000.;
}

double Constant::get_spectral_density(double wavelength) {
    return this->value;
}

Constant::Constant(double value) {
    this->value = value;
}

Constant::Constant() {
    this->value = 1.0;

}

IdealEtalon::IdealEtalon(double d, double n, double theta, double R) : d(d/1000.), n(n), theta(theta), R(R) {
    this->cF = this->coefficient_of_finesse(R);
}

double IdealEtalon::coefficient_of_finesse(double R) {
    return 4.*R / ((1.-R)*(1.-R));
}

double IdealEtalon::FSR() {
    return 299792458. /  (2.*n*d*cos(theta));
}

double IdealEtalon::T(double wl, double theta, double d, double n, double cF) {
    //delta = (2. * math.pi / wl) * 2. * n * math.cos(theta) * d
    //return 1. / (1. + cF * np.sin(0.5 * delta) ** 2)

    double delta = (2. * M_PI * n * cos(theta) * d) / wl ;
    double sind = sin(delta);
    return 1. / (1. + cF * sind*sind);
}

double IdealEtalon::get_spectral_density(double wavelength) {
    double res = this->T(wavelength/1E6, theta, d, n, cF);
    return res;
}


Blackbody::Blackbody(double T): T(T){};

double Blackbody::planck(const double& T, const double& wavelength) {
    const double 	hPlanck = 6.62606896e-34;
    const double 	speedOfLight = 2.99792458e8;
    const double 	kBoltzmann = 1.3806504e-23;
    double a = 2.0 * hPlanck * speedOfLight*speedOfLight;
    double b = hPlanck * speedOfLight / (wavelength * kBoltzmann * T);
    double intensity = a / (pow(wavelength,5) * (exp(b) - 1.0));
    return intensity;
}

double Blackbody::get_spectral_density(double wavelength) {
    return this->planck(this->T, wavelength/1E6);
}


PhoenixSpectrum::PhoenixSpectrum(std::string spectrum_file, std::string wavelength_file, const double &min_wavelength,
                                 const double &max_wavelength) {
    this->read_spectrum(spectrum_file, wavelength_file, min_wavelength, max_wavelength);
}



void PhoenixSpectrum::read_spectrum(std::string spectrum_file, std::string wavelength_file,
                                    const double &min_wavelength, const double &max_wavelength) {
// read in wavelength file
    std::unique_ptr<CCfits::FITS> ptr_FITS_wl(new CCfits::FITS(wavelength_file, CCfits::Read, true));
    CCfits::PHDU& wl = ptr_FITS_wl->pHDU();
    std::valarray<double>  contents_wl;
    wl.readAllKeys();
    wl.read(contents_wl);

    // this doesn't print the data, just header info.
    std::cout << wl << std::endl;

    long ax1(wl.axis(0));
    long min_idx=ax1;
    long max_idx=0;
    double min_wavelength_angstrom = min_wavelength*10000.;
    double max_wavelength_angstrom = max_wavelength*10000.;

    // find wavelength limits
    for (long j = 0; j < ax1; ++j)
    {
        if (contents_wl[j] < min_wavelength_angstrom)
            min_idx = j;
        if ( contents_wl[j] > max_wavelength_angstrom)
            max_idx = j;
    }
    std::cout << min_idx << "\t" << max_idx ;

    std::unique_ptr<CCfits::FITS> ptr_FITS_spectrum(new CCfits::FITS(spectrum_file, CCfits::Read, true));
    CCfits::PHDU& spec = ptr_FITS_spectrum->pHDU();
    std::valarray<double>  contents_spec;
    spec.readAllKeys();
    spec.read(contents_spec);
    // find wavelength limits
    for (long j = min_idx; j < max_idx; ++j) {
        this->data[contents_wl[j]/10000.] = contents_spec[j]/1E14;
    }
}

double PhoenixSpectrum::get_spectral_density(double wavelength) {
    return interpolate(this->data, wavelength);
    /*
    long idx = 0;
    while (idx < this->wavelength.size() && this->wavelength[idx] < wavelength) {
        ++idx;
    }
    double frac = (wavelength - this->wavelength[idx-1]) / (this->wavelength[idx] - this->wavelength[idx-1]);
    return lerp(this->spectrum[idx-1], this->spectrum[idx], frac);
     */
}

LineList::LineList(std::string linelist) {
    this->read_spectrum(linelist);

}



void LineList::read_spectrum(std::string linelist) {
    std::ifstream       file(linelist.c_str());

    for(CSVIterator loop(file); loop != CSVIterator(); ++loop)
    {
        this->data.insert(std::pair<double, double> (std::stod((*loop)[0]), std::stod((*loop)[1])));
    }

}

std::vector<double> LineList::get_spectrum(std::vector<double> wavelength){
    std::vector<double> spectrum;
    for(auto const & wl: wavelength){
        double d = fabs(wl - this->data.begin()->first);
        for(auto iterator = ++this->data.begin(); iterator!= this->data.end(); ++iterator)
        {
            if (fabs(wl - iterator->first) < d){
                d = wl - iterator->first;
            }
            else
            {
                spectrum.push_back(iterator->second);
                break;
            }
        }
    }
    return spectrum;
}

double LineList::get_spectral_density(double wavelength) {

}

std::vector<double> LineList::get_wavelength() {
    std::vector<double> wavelength;
    for (auto m: this->data){
        wavelength.push_back(m.first);
    }
    return wavelength;
}