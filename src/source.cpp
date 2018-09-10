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
//#include <Histogram.h>

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
    return std::vector<double>(wavelength.size());
}

double Source::get_spectral_density(double wavelength) {
    return 1.;
}
std::vector<float> Source::get_spectrum(std::vector<double> wavelength) {
    if(this->mode) {
        std::vector<float> spectrum;
        std::vector<double> diff;

        for (int i=0; i < wavelength.size(); i++) {
            diff.push_back(this->shift * (wavelength[i + 1] - wavelength[i]));
        }
        diff.push_back(diff.back());

        for (int i=0; i < wavelength.size(); i++) {
            double a = this->shift * wavelength[i] - diff[i] / 2.;
            double b = this->shift * wavelength[i] + diff[i] / 2.;
            float result = float(this->integral_s(a, b, this->integration_steps) / diff[i]);
            spectrum.push_back(result);
        }

        return spectrum;
    }
    else{
        std::vector<float> spectrum;

        for (int i = 0; i < wavelength.size(); i++) {
            std::cout<<spectrum[i]<<std::endl;
            spectrum.push_back(float(get_spectral_density(wavelength[i])));
        }

        return spectrum;
    }
}

void Source::set_integration_steps(int n) {
    this->integration_steps = n;
}

double Source::integral_s(double a, double b, int n) {
        double step = (b-a) / n;  // width of each small rectangle
        double area = 0.0;  // signed area
        for (int i = 0; i < n; i ++) {
            area += this->get_spectral_density(a + (i + 0.5) * step) * step; // sum up each small rectangle
        }
        return area;
}

void Source::set_doppler_shift(double shift) {
    this->shift = 1. + shift / 300000000.;
}

void Source::scale_spectral_density() {

    double a = min_w;
    double b = max_w;
    double inten_pho = 5.03E10;
    int n =1000000;
    double step = (b-a) / n;  // width of each small rectangle
    double area = 0.0;  // signed area
    for (int i = 0; i < n; i ++) {
        area += this->get_spectral_density(a + (i + 0.5) * step) * (inten_pho) * (a + (i + 0.5) * step) * (step); // sum up each small rectangle
    }

    s_val = s_val * pow(10, mag/(-2.5))*v_zp / (area);
    //std::cout<<s_val<<std::endl;

}

double Constant::get_spectral_density(double wavelength) {
    return this->value;
}

Constant::Constant(double value, double min_w, double max_w) {
    this->value = value; // uW per um (micro watts per micro meter// )
    this -> min_w = min_w;
    this -> max_w = max_w;
}

Constant::Constant() {
    this->value = 0.01 ; // uW per um (micro watts per micro meter)
    min_w = 0;
    max_w = 1000000;
}

Blackbody::Blackbody(double T, double mag): T(T){

    this -> mag = mag;

    min_w = 0; // Maybe set a cutoff based off intensity?
    max_w = 1000000;

    scale_spectral_density();

}

double Blackbody::planck(const double& T, const double& wavelength) {
    const double 	hPlanck = 6.62606896e-34; //J * s;
    const double 	speedOfLight = 2.99792458e8; //m / s;
    const double 	kBoltzmann = 1.3806504e-23; //J / K
    double a = 2.0 * hPlanck * speedOfLight*speedOfLight;
    double b = hPlanck * speedOfLight / (wavelength * kBoltzmann * T);
    double c = 1E0; // Conversion factor for intensity
    double intensity = (s_val) * (c) * a / (pow(wavelength,5) * (exp(b) - 1.0)); // (J / s ) / m^3 -> (uJ / s) / ( m^2 * um )
    return intensity;
}

double Blackbody::get_spectral_density(double wavelength) {
    return this->planck(this->T, wavelength/1E6);
}

PhoenixSpectrum::PhoenixSpectrum(std::string spectrum_file, std::string wavelength_file, double mag) {
    this->read_spectrum(spectrum_file, wavelength_file);
    this->mag = mag;
    scale_spectral_density();
}

void PhoenixSpectrum::read_spectrum(std::string spectrum_file, std::string wavelength_file) {
// read in wavelength file
    std::unique_ptr<CCfits::FITS> ptr_FITS_wl(new CCfits::FITS(wavelength_file, CCfits::Read, true));
    CCfits::PHDU& wl = ptr_FITS_wl->pHDU();
    std::valarray<double>  contents_wl;
    wl.readAllKeys();
    wl.read(contents_wl);

    long ax1(wl.axis(0));
    long min_idx=ax1;
    long max_idx=0;
    double min_wavelength_angstrom = min_w*10000.;
    //std::cout<<min_wavelength;
    double max_wavelength_angstrom = max_w*10000.;

    // find wavelength limits
    for (long j = 0; j < ax1; ++j)
    {
        if (contents_wl[j] < min_wavelength_angstrom)
            min_idx = j;
        if ( contents_wl[j] > max_wavelength_angstrom)
            max_idx = j;
    }
    std::cout << "min_idx"<<":"<<min_idx<<"\t";
    std::cout<<"max_idx"<<":"<<max_idx<<"\t";

    std::unique_ptr<CCfits::FITS> ptr_FITS_spectrum(new CCfits::FITS(spectrum_file, CCfits::Read, true));
    CCfits::PHDU& spec = ptr_FITS_spectrum->pHDU();
    std::valarray<double>  contents_spec;
    spec.readAllKeys();
    spec.read(contents_spec);
    // find wavelength limits
    // convert contents_spec from erg/s/cm^2/cm to uW/s/m^2/um
    for (long j = min_idx; j < max_idx; ++j) {
        this->data[contents_wl[j]/10000.] = (1E3) * contents_spec[j];
        //std::cout<<contents_spec[j]<<":";
    }

}

double PhoenixSpectrum::get_spectral_density(double wavelength) {

    return (s_val)*interpolate(this->data, wavelength);
    //std::cout<<a<<":";

    /*
    long idx = 0;
    while (idx < this->wavelength.size() && this->wavelength[idx] < wavelength) {
        ++idx;
    }
    double frac = (wavelength - this->wavelength[idx-1]) / (this->wavelength[idx] - this->wavelength[idx-1]);
    return lerp(this->spectrum[idx-1], this->spectrum[idx], frac);
     */
}

CoehloSpectrum::CoehloSpectrum(std::string spectrum_file, double mag) {
    this->read_spectrum(std::move(spectrum_file));
    this->mag = mag;
    scale_spectral_density();
}

void CoehloSpectrum::read_spectrum(std::string spectrum_file) {
//// read in wavelength file
//    std::unique_ptr<CCfits::FITS> ptr_FITS_wl(new CCfits::FITS(wavelength_file, CCfits::Read, true));
//    CCfits::PHDU& wl = ptr_FITS_wl->pHDU();
//    std::valarray<double>  contents_wl;
//    wl.readAllKeys();
//    wl.read(contents_wl);
//
//    long ax1(wl.axis(0));
//    long min_idx=ax1;
//    long max_idx=0;
//    double min_wavelength_angstrom = min_w*10000.;
//    //std::cout<<min_wavelength;
//    double max_wavelength_angstrom = max_w*10000.;
//
//    // find wavelength limits
//    for (long j = 0; j < ax1; ++j)
//    {
//        if (contents_wl[j] < min_wavelength_angstrom)
//            min_idx = j;
//        if ( contents_wl[j] > max_wavelength_angstrom)
//            max_idx = j;
//    }
//    std::cout << "min_idx"<<":"<<min_idx<<"\t";
//    std::cout<<"max_idx"<<":"<<max_idx<<"\t";

    std::unique_ptr<CCfits::FITS> ptr_FITS_spectrum(new CCfits::FITS(spectrum_file, CCfits::Read, true));
    CCfits::PHDU& spec = ptr_FITS_spectrum->pHDU();
    std::valarray<double>  contents_spec;
    spec.readAllKeys();
    spec.read(contents_spec);
    long ax1(spec.axis(0));
    long max_idx=ax1;
    long min_idx=0;
    std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(spectrum_file, CCfits::Read));
    double minwl;
    pInfile->pHDU().readKey("CRVAL1", minwl);

    double wld;
    pInfile->pHDU().readKey("CDELT1", wld);

    // find wavelength limits
    // convert contents_spec from erg/s/cm^2/A to uW/s/m^2/um
    for (long j = min_idx; j < max_idx; ++j) {
        this->data[(minwl+j*wld)/10000.] = (1E11) * contents_spec[j];
        //std::cout<<contents_spec[j]<<":";
    }

}

double CoehloSpectrum::get_spectral_density(double wavelength) {

    return (s_val)*interpolate(this->data, wavelength);

}

CustomSpectrum::CustomSpectrum(const std::string spectrum_file, double min_w, double max_w, double mag){

    this -> min_w = min_w;
    this -> max_w = max_w;
    this -> mag = mag;

    read_spectrum(spectrum_file);
    scale_spectral_density();

}

CustomSpectrum::CustomSpectrum(const std::string spectrum_file, std::string wave_file, double mag){

    this -> mag = mag;

    read_spectrum(spectrum_file, wave_file);
    scale_spectral_density();

}

void CustomSpectrum::read_spectrum(const std::string spectrum_file) {
//// read in wavelength file

    std::unique_ptr<CCfits::FITS> ptr_FITS_spectrum(new CCfits::FITS(spectrum_file, CCfits::Read, true));
    CCfits::PHDU& spec = ptr_FITS_spectrum->pHDU();
    std::valarray<double>  contents_spec;
    spec.readAllKeys();
    spec.read(contents_spec);
    long ax1(spec.axis(0));
    long max_idx=ax1;
    long min_idx=0;
    std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(spectrum_file, CCfits::Read));

    // find wavelength limits
    // convert contents_spec from erg/s/cm^2/A to uW/s/m^2/um
    for (long j = min_idx; j < max_idx; ++j) {
        this->data[( (max_w-min_w)*j / max_idx + min_w ) /10000.] = (1E11) * contents_spec[j];
        //std::cout<<contents_spec[j]<<":";
    }

}

void CustomSpectrum::read_spectrum(const std::string spectrum_file, std::string wave_file) {
//// read in wavelength file
    std::unique_ptr<CCfits::FITS> ptr_FITS_wl(new CCfits::FITS(wave_file, CCfits::Read, true));
    CCfits::PHDU& wl = ptr_FITS_wl->pHDU();
    std::valarray<double>  contents_wl;
    wl.readAllKeys();
    wl.read(contents_wl);

    long ax1(wl.axis(0));
    long max_idx=ax1;
    long min_idx=0;
    double min_wavelength_angstrom = min_w*10000.;
    //std::cout<<min_wavelength;
    double max_wavelength_angstrom = max_w*10000.;

    // find wavelength limits
    for (long j = 0; j < ax1; ++j)
    {
        if (contents_wl[j] < min_wavelength_angstrom)
            min_idx = j;
        if ( contents_wl[j] > max_wavelength_angstrom)
            max_idx = j;
    }
    std::cout << "min_idx"<<":"<<min_idx<<"\t";
    std::cout<<"max_idx"<<":"<<max_idx<<"\t";

    std::unique_ptr<CCfits::FITS> ptr_FITS_spectrum(new CCfits::FITS(spectrum_file, CCfits::Read, true));
    CCfits::PHDU& spec = ptr_FITS_spectrum->pHDU();
    std::valarray<double>  contents_spec;
    spec.readAllKeys();
    spec.read(contents_spec);
    std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(spectrum_file, CCfits::Read));

    // find wavelength limits
    // convert contents_spec from erg/s/cm^2/A to uW/s/m^2/um
    for (long j = min_idx; j < max_idx; ++j) {
        this->data[ contents_wl[j] / 10000.] = (1E11) * contents_spec[j];
        //std::cout<<contents_spec[j]<<":";
    }

}

double CustomSpectrum::get_spectral_density(double wavelength) {

    return (s_val)*interpolate(this->data, wavelength);

}

LineList::LineList(const std::string linelist_file, double scaling) {
    mode = false;
    this -> scaling = scaling;
    this->read_spectrum(linelist_file);

    smooth_spectrum();

}

void LineList::read_spectrum(std::string linelist_file) {

//    std::ifstream       file(linelist_file.c_str());

//    for(CSVIterator loop(file); loop != CSVIterator(); ++loop)
//    {
//        std::cout<<"Test"<<std::endl;
//        event.push_back(stod((*loop)[0]));
//        intensity.push_back(stod((*loop)[1]));
//        this->data.insert(std::pair<double, double> (stod((*loop)[0]), stod((*loop)[1])));
//
//    }

    std::ifstream open(linelist_file.c_str(), std::ios::in);

    if (!open) {
        std::cout << "No File Found";
    }

    double a, c;
    char b;
    int i = 0;

    while (open >> a >> b >> c) {
        i = i + 1;

        event.push_back(a);

        intensity.push_back(scaling * c);

        this->data.insert(std::pair<double, double> (a, scaling * c));

    }

    length = i;



}

void LineList::smooth_spectrum(){



}

std::vector<float> LineList::get_spectrum(std::vector<double> wavelength){
    std::vector<float> spectrum;
    for(auto const & wl: wavelength){

        spectrum.push_back(float(data[wl]));

    }
    return spectrum;
}

double LineList::get_spectral_density(double wavelength) {

    return data[wavelength];

}

std::vector<double> LineList::get_wavelength() {
    std::vector<double> wavelength;
    for (auto m: this->data){
        wavelength.push_back(m.first);
    }
    return wavelength;
}

template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}