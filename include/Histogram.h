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

/*!
 * \class Histogram
 * \brief The general statistical object that spectra is a child of
 *
 * This class allows for a variety of statistical computations to be carried out
 * on members from source. It represents sources as a histogram with a cdf and sampling
 * methods.
 *
 */

class Histogram{
protected:

    void Read_Data(string file_path);
    void Read_Data(string path_1, string path_2);
    void Create_cdf();

public:
    unsigned long length=0; //This is the number of elements in the histogram object
    int resolution=0;
    /// Cumalative distribution function for the histogram
    vector<double> cdf;
    /// The event space of the histogram (wavelengths)
    vector<double> event;
    /// The number of events
    vector<double> intensity;

    /// Differential of event
    vector<double> d_event;
    /// Differential of intensity
    vector<double> d_intensity;
    /// Differential of the cdf (the pdf)
    vector<double> d_cdf; //This is the same as a pdf vector

    Histogram();
    /// Generates a histogram from a csv
    Histogram(string file_path); //Generate histogram from file path
    /// Generates a histogram from event csv and intensity csv
    Histogram(string path_1, string path_2);

    /// Generates Histogram from a vector of events and weights
    Histogram(vector<double> events, vector<double> weights); //Generate histogram manually
    /// Will resample another histogram with a number of elements
    Histogram(int num,Histogram &histogram); //Generate sample histogram from original with custom number of elements
    /// Resamples a histogram with a different resolution and number of elements
    Histogram(int num,int res, Histogram &histogram);//Same as above, but also allows for res x the sampling resolution of the original histogram

    /*!
     *Returns a test statistic regarding whether or not two histograms are similar
     *Null hypothesis: The histograms are drawn from the same source.
     *We want the test statistic to be less than c(a)*sqrt((n+m)/(n*m)) to keep the null hypothesis
     *...Otherwise we reject it
     *a|    0.10	0.05	0.025	0.01	0.005	0.001
     *c(a)| 1.22	1.36	1.48	1.63	1.73	1.95
     *
     * @param sim_histogram the comparison histogram
     * @return
     */
    double Kolmogorov_Smirnov_test(Histogram &sim_histogram);

    /*!
    * Allows wavelengths to be sampled from a histogram using
    * CDF sampling methods
    *
    * @param sample_value  A value in the range 0 to 1
    */
    template<class T> double Sample(T sample_value); //Sample using point methods

    /// Allows for multi-sampling
    template<class T> vector<double> Sample(vector<T> sample_value);

    /// Allows for sampling using a linear correction to the CDF
    template<class T> double Linear_sample(T sample_value); //Sample using point methods and then correct error to first order

    /// Allows for multi-sampling with linear correction
    template<class T> vector<double> Linear_sample(vector<T> sample_value);

    bool mode = true;

};

Histogram::Histogram() {
}

Histogram::Histogram(string file_path){
    Read_Data(file_path);
    Create_cdf();
}

Histogram::Histogram(string path_1, string path_2){
    Read_Data(path_1,path_2);
    Create_cdf();
}

Histogram::Histogram(vector<double> events, vector<double> weights){

    event=events;
    intensity=weights;
    length=events.size();
    Create_cdf();

}

//This constructor is different from the rest
//It takes another object in the class and then creates a new object of a different size
Histogram::Histogram(int num,Histogram &histogram){

    vector<double> emp_weights(histogram.length);

    long long int seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed);
    std::uniform_real_distribution<double> uniform_distribution (0.0,1.0);

    for(int i=0; i<num; i++) {
        double sample_value = uniform_distribution(generator);
        long dist=distance(histogram.cdf.begin(),lower_bound(histogram.cdf.begin(),histogram.cdf.end(),sample_value))-1;
        if(dist == -1){dist=0;};

        emp_weights[dist]=emp_weights[dist]+1;

        //cout<<spectra.Linear_sample(sample_value)<<":"<<sample_value<<endl;
    }

    event=histogram.event;
    intensity=emp_weights;
    length=histogram.length;
    Create_cdf();

}

Histogram::Histogram(int num,int res, Histogram &histogram){
    //more resolution means more elements
    length=histogram.length*(res);
    vector<double> emp_weights(static_cast<unsigned>(histogram.length*res));
    vector<double> emp_lambdas(static_cast<unsigned>(histogram.length*res));

    for(int i=0; i<length; i++){
        emp_weights[i]=0;
        //We have to resize the event space i.e the lambdas
        emp_lambdas[i]=histogram.event[0]+static_cast<double>(i-1) * (histogram.event[histogram.length-1]-histogram.event[0]) / static_cast<double>(length);
    }

    int place;
    resolution=res;

    long long int seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed);
    std::uniform_real_distribution<double> uniform_distribution (0.0,1.0);

    for(int i=0; i<num; i++) {
        double sample_value = uniform_distribution(generator);
        //CDF[0] is not meaningful so we subtract 1 to discard it with an iff statement
        long dist = distance(histogram.cdf.begin(), lower_bound(histogram.cdf.begin(), histogram.cdf.end(), sample_value))-1;
        if(dist == -1){dist=0;};

        //finds the new index for wheights given the new resolution
        place=static_cast<int>((dist)*res+floor(res*(histogram.cdf[dist+1]-sample_value)/(histogram.cdf[dist+1]-histogram.cdf[dist])));
        emp_weights[place]=emp_weights[place]+1;

        //cout<<emp_weights[place]<<":"<<emp_lambdas[place]<<endl;

        //cout<<spectra.Linear_sample(sample_value)<<":"<<sample_value<<endl;
    }

    event=emp_lambdas;
    intensity=emp_weights;
    Create_cdf();

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

void Histogram::Create_cdf(){

        cdf.push_back(0);

        for (int i = 0; i < length - 1; i++) {
            d_event.push_back(event[i + 1] - event[i]);
            d_intensity.push_back(intensity[i + 1] - intensity[i]);
            //The factor of 1/2 comes from the fact we are numerically integrating and thus
            //taking the average of intensity at two nearby points
            cdf.push_back(cdf[i] + (0.5) * (d_event[i]) * (intensity[i + 1] + intensity[i]));
            d_cdf.push_back((0.5) * (d_event[i]) * (intensity[i + 1] + intensity[i]));
        }

        d_event.push_back(event[length - 1] - event[length - 2]);
        d_cdf.push_back((0.5) * (d_event[length - 1]) * (intensity[length - 1] + intensity[length - 2]));
        d_intensity.push_back(intensity[length - 1] - intensity[length - 2]);

        d_event.push_back(event[length - 1] - event[length - 2]);
        d_cdf.push_back((0.5) * (d_event[length - 1]) * (intensity[length - 1] + intensity[length - 2]));
        d_intensity.push_back(intensity[length - 1] - intensity[length - 2]);

        double norm = cdf[length - 1];
        //We have to create a normalized cdf
        //So we have to go back through and divide by the integration constant
        //This is expensive...however note that cdf is already ordered for us
        for (int i = 0; i < length; i++) {
            cdf[i] = cdf[i] / norm;
            d_cdf[i] = d_cdf[i] / norm;
        }

    return;
}

//There's going to be an issue if the the event elements of the two objects aren't the same.
//This test works by finding the sup of the difference between the two histograms cdfs
// then comparing it to a test statistic
double Histogram::Kolmogorov_Smirnov_test(Histogram &sim_histogram){

    double max_diff=0;
    double diff=0;
    long dist_1;
    long dist_2;
    for(double k=event[0]; k<event[length-1]; k=k+(event[length-1]-event[0])/max(length,sim_histogram.length)){

        //we need two corresponding distances since the histograms may have different numbers of elements
        dist_1=distance(event.begin(),lower_bound(event.begin(), event.end(), k));
        dist_2=distance(sim_histogram.event.begin(),lower_bound(sim_histogram.event.begin(), sim_histogram.event.end(), k));

        //k may run outside the histogram range of lambdas
        if(dist_1>length-1){dist_1=length-1;};
        if(dist_2>sim_histogram.length-1){dist_2=sim_histogram.length-1;};

        diff=cdf[dist_1]-sim_histogram.cdf[dist_2];

        if(abs(diff)>max_diff){max_diff=abs(diff);};

    }

    return max_diff;
}

template<class T>
double Histogram::Sample(T sample_value){

        long dist = distance(cdf.begin(), lower_bound(cdf.begin(), cdf.end(), sample_value)) - 1;
        if (dist == -1) { dist = 0; };
        //There's a lot going on here. distance() finds the number of elements seperating two values of cdf
        //lower_bound finds the closest member of cdf to sample_valuepl

        return event[dist];

}


template<class T>
vector<double> Histogram::Sample(vector<T> sample_value){
    vector<double> vec;
    for(int i=0; i<sample_value.size(); i++){vec.push_back(Sample(sample_value[i]));}
    return vec;
}



template<class T>
double Histogram::Linear_sample(T sample_value) {
    long dist = distance(cdf.begin(), lower_bound(cdf.begin(), cdf.end(), sample_value))-1;
    if(dist == -1){dist=0;};

    return event[dist] + (d_event[dist]) / (d_cdf[dist+1]) * (cdf[dist+1] - sample_value);
}

template<class T>
vector<double> Histogram::Linear_sample(vector<T> sample_value){
    vector<double> vec;
    for(int i=0; i<sample_value.size(); i++){vec.push_back(Linear_sample(sample_value[i]));}
    return vec;
}

#endif //_HISTOGRAM_H
