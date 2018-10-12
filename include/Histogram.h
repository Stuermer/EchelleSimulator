//
// Created by zacha_000 on 6/14/2017.
//

#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <iostream>
#include <vector>
#include <chrono>
#include <vector>
#include <random>
#include <fstream>
#include "random_generator.h"

using namespace std;

class Histogram{
protected:

    void Read_Data(string file_path);
    void Read_Data(string path_1, string path_2);

public:
    unsigned long length=0; //This is the number of elements in the histogram object
    int resolution=0;
    vector<double> cdf;
    vector<double> event;
    vector<double> intensity;
    vector<double> flux;

    vector<double> d_flux;
    vector<double> d_event;
    vector<double> d_intensity;
    vector<double> d_cdf; //This is the same as a pdf vector

    //std::piecewise_constant_distribution<float> distribution;

    Histogram();
    Histogram(string file_path); //Generate histogram from file path
    Histogram(string path_1, string path_2);
    virtual double Calc_flux() {return 1.;};

    Histogram(vector<double> events, vector<double> weights); //Generate histogram manually

    bool mode = true;

    double v_zp=8660006000.0; //The reference flux is obtained by integrating vega
    // over a bessel filter and has units photons/m^2/s

};

Histogram::Histogram() {
}

Histogram::Histogram(string file_path){
    Read_Data(file_path);
}

Histogram::Histogram(string path_1, string path_2){
    Read_Data(path_1,path_2);
}

Histogram::Histogram(vector<double> events, vector<double> weights){

    event=events;
    intensity=weights;
    length=events.size();

}

void Histogram::Read_Data(string file_path) {

    ifstream open(file_path, ios::in);

    if (!open) {
        cout << "No File Found";
    }

    double a, b;
    int i = 0;

    while (open >> a >> b) {
        i = i + 1;

        event.push_back(a);

        intensity.push_back(b);

    }

    length = i;

}

void Histogram::Read_Data(string path_1, string path_2) {

    ifstream open_1(path_1, ios::in);

    if (!open_1) {
        cout << "No File Found";
    }

    double a, b;
    int i = 0;
    int j = 0;

    while (open_1 >> a) {
        i = i + 1;

        event.push_back(a);

    }

    length = i;

    ifstream open_2(path_2, ios::in);

    if (!open_2) {
        cout << "No File Found";
    }

    while (open_2 >> b) {
        j = j + 1;

        intensity.push_back(b);

    }

    if(length!=j){ cout<<"Data Mismatch"<<endl; }

    return;

}

class Smooth_Spectra: public Histogram {
public:
    Smooth_Spectra() {
    }

    Smooth_Spectra(vector<double> events, vector<double> weights) : Histogram(events,weights) {
        Create_flux();
        Create_cdf();
    }

    void Create_cdf();
    void Create_flux();
    double Calc_flux();

    template<class T>
    double Sample(T sample_value); //Sample using point methods

    template<class T>
    vector<double> Sample(vector<T> sample_value);

    template<class T>
    double Linear_sample(T sample_value); //Sample using point methods and then correct error to first order

    template<class T>
    vector<double> Linear_sample(vector<T> sample_value);
};

void Smooth_Spectra::Create_cdf(){

    cdf.push_back(0);

    for (int i = 0; i < length - 1; i++) {
        d_event.push_back(event[i + 1] - event[i]);
        d_flux.push_back(flux[i + 1] - flux[i]);
        //The factor of 1/2 comes from the fact we are numerically integrating and thus
        //taking the average of intensity at two nearby points
        cdf.push_back(cdf[i] + (0.5) * (d_event[i]) * (flux[i + 1] + flux[i]));
        d_cdf.push_back((0.5) * (d_event[i]) * (flux[i + 1] + flux[i]));
    }

    d_event.push_back(event[length - 1] - event[length - 2]);
    d_cdf.push_back((0.5) * (d_event[length - 1]) * (flux[length - 1] + flux[length - 2]));
    d_intensity.push_back(flux[length - 1] - flux[length - 2]);

    d_event.push_back(event[length - 1] - event[length - 2]);
    d_cdf.push_back((0.5) * (d_event[length - 1]) * (flux[length - 1] + flux[length - 2]));
    d_intensity.push_back(flux[length - 1] - flux[length - 2]);

    double norm = cdf[length - 1];
    //We have to create a normalized cdf
    //So we have to go back through and divide by the integration constant
    //This is expensive...however note that cdf is already ordered for us
    for (int i = 0; i < length; i++) {
        cdf[i] = cdf[i] / norm;
        d_cdf[i] = d_cdf[i] / norm;
    }

    //std::piecewise_constant_distribution<float> distribution (event.begin(),event.end(),intensity.begin());

    return;
}

void Smooth_Spectra::Create_flux(){

    //0.503 is for assuming intensity is erg/s/cm^2/cm
    //and that event~(wavelength) is in A

    //5.03*10^10 is for assuming intensity is uW/m^2/um
    //and that event is in um

    double ch_factor = 5.03E10;

    for (int i = 0; i < length; i++) {
        flux.push_back(intensity[i] * (ch_factor) * (event[i]));

    }

    flux.push_back(intensity[length - 1] * (ch_factor) * event[length - 1]);

    return;
}

double Smooth_Spectra::Calc_flux(){

    double tflux = 0;

    for (int i = 0; i < length; i++) {

        tflux = tflux + (0.5) * (flux[i + 1] + flux[i]) * (d_event[i]);

    }

    return tflux;

}


template<class T>
double Smooth_Spectra::Sample(T sample_value){

    long dist = distance(cdf.begin(), lower_bound(cdf.begin(), cdf.end(), sample_value)) - 1;
    if (dist == -1) { dist = 0; };
    //There's a lot going on here. distance() finds the number of elements seperating two values of cdf
    //lower_bound finds the closest member of cdf to sample_valuepl

    return event[dist];

}

template<class T>
vector<double> Smooth_Spectra::Sample(vector<T> sample_value){
    vector<double> vec;
    for(int i=0; i<sample_value.size(); i++){vec.push_back(Sample(sample_value[i]));}
    return vec;
}



template<class T>
double Smooth_Spectra::Linear_sample(T sample_value) {
    long dist = distance(cdf.begin(), lower_bound(cdf.begin(), cdf.end(), sample_value))-1;
    if(dist == -1){dist=0;};

    return event[dist] + (d_event[dist]) / (d_cdf[dist+1]) * (cdf[dist+1] - sample_value);
}

template<class T>
vector<double> Smooth_Spectra::Linear_sample(vector<T> sample_value){
    vector<double> vec;
    for(int i=0; i<sample_value.size(); i++){vec.push_back(Linear_sample(sample_value[i]));}
    return vec;
}

class Discrete_Spectra: public Histogram {
public:
    Discrete_Spectra() {

    }

    Discrete_Spectra(vector<double> events, vector<double> weights) : Histogram(events,weights) {
        Create_flux();
        Create_cdf();
        mode = 0;
    }

    void Create_cdf();
    void Create_flux();
    double Calc_flux();

    template<class T>
    double Sample(T sample_value); //Sample using point methods

    template<class T>
    vector<double> Sample(vector<T> sample_value);

};

void Discrete_Spectra::Create_cdf(){

    cdf.push_back(0);

    for (int i = 0; i < length - 1; i++) {
        d_flux.push_back(flux[i + 1] - flux[i]);
        //The factor of 1/2 comes from the fact we are numerically integrating and thus
        //taking the average of intensity at two nearby points
        cdf.push_back(cdf[i] + flux[i]);
        d_cdf.push_back(flux[i]);
    }

    d_cdf.push_back(flux[length - 2]);
    d_intensity.push_back(flux[length - 2]);

    d_cdf.push_back(flux[length - 1]);
    d_intensity.push_back(flux[length - 1]);

    double norm = cdf[length - 1];
    //We have to create a normalized cdf
    //So we have to go back through and divide by the integration constant
    //This is expensive...however note that cdf is already ordered for us
    for (int i = 0; i < length; i++) {
        cdf[i] = cdf[i] / norm;
        d_cdf[i] = d_cdf[i] / norm;
    }

    //std::piecewise_constant_distribution<float> distribution (event.begin(),event.end(),intensity.begin());

    return;
}

void Discrete_Spectra::Create_flux(){

    //0.503 is for assuming intensity is erg/s/cm^2/cm
    //and that event~(wavelength) is in A

    //5.03*10^10 is for assuming intensity is uW/m^2/um
    //and that event is in um

    // This needs to be one sense the spectral lines are thin.
    double ch_factor = 1.;

    for (int i = 0; i < length; i++) {
        flux.push_back(intensity[i] * (ch_factor) * (event[i]));
    }

    flux.push_back(intensity[length - 1] * (ch_factor) * event[length - 1]);

    return;
}

double Discrete_Spectra::Calc_flux(){

    double tflux = 0;


    for (int i = 0; i < length; i++) {
        tflux = tflux + flux[i];

    }

    return tflux;

}

template<class T>
double Discrete_Spectra::Sample(T sample_value){

    long dist = distance(cdf.begin(), lower_bound(cdf.begin(), cdf.end(), sample_value)) - 1;
    if (dist == -1) { dist = 0; };
    //There's a lot going on here. distance() finds the number of elements seperating two values of cdf
    //lower_bound finds the closest member of cdf to sample_valuepl

    return event[dist];

}

template<class T>
vector<double> Discrete_Spectra::Sample(vector<T> sample_value){
    vector<double> vec;
    for(int i=0; i<sample_value.size(); i++){vec.push_back(Sample(sample_value[i]));}
    return vec;
}

#endif //_HISTOGRAM_H
