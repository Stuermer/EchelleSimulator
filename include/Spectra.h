//
// Created by zrobertson on 6/27/17.
//

#ifndef RANDOMGENERATOR_SPECTRA_H
#define RANDOMGENERATOR_SPECTRA_H

#include "Histogram.h"
#include <iostream>

class Spectra:public Histogram{

public:
    Spectra();
    Spectra(string path);
    Spectra(string path_1, string path_2);

    Spectra(vector<double> m_event, vector<double> m_intensity);
    Spectra(Spectra spectra, Histogram filter);
    //Filters the spectra and adds the result as an object
    //Technically deprecated because of Pass_filter

    const double v_zp=866000.0; //The reference flux is obtained by integrating vega
                                // over a bessel filter and has units photons/cm^2/s

    double magnitude;

    vector<double> dflux; //Units are flux per wavelength  photons/s/cm^2/A

    void Create_dflux();
    void Set_magnitude(double mag, double t_flux);

    Spectra Pass_filter(Histogram filter); //Passes the filter over Spectra
    vector<double> Filter_dflux(Histogram filter);

    double Relative_flux(Histogram m_order, vector<Histogram> m_orders, double r_flux);
    //Calculates flux through m_order assuming
    // there are m_orders and that r_flux is conserved (slower by a factor of num(m_orders))

    double Relative_flux(Histogram m_order, vector<Histogram> m_orders); //Same as above but assumes r_flux = total_flux

    double Calc_flux(const Histogram& filter);
    double Calc_flux();

    double Photon_count(double m_time, double m_area);
    double Exp_time(double s_noise);

};

Spectra::Spectra(){

}

Spectra::Spectra(string path):Histogram(path){

    Create_dflux();
    Calc_flux();
}

Spectra::Spectra(string path_1, string path_2):Histogram(path_1, path_2){

    Create_dflux();
    Calc_flux();
}

Spectra::Spectra(vector<double> m_event, vector<double> m_intensity):Histogram(m_event, m_intensity){

    Create_dflux();
    Calc_flux();
}

Spectra::Spectra(Spectra spectra, Histogram filter):Histogram(){

    long dist;
    double weight;
    vector<double> t_intensity;

    for(int i=0; i<spectra.length; i++){

        dist=distance(filter.event.begin(),lower_bound(filter.event.begin(),filter.event.end(),spectra.event[i]));
        weight=filter.intensity[dist]+(filter.d_intensity[dist]) / (filter.d_event[dist]) * (filter.event[dist] - spectra.event[i]);

        if(filter.event[spectra.length] - spectra.event[i]<0){weight=0;}
        if(filter.event[0] - spectra.event[i]>0){weight=0;}

        t_intensity.push_back(spectra.intensity[i]*weight);

    }

    event=spectra.event;
    intensity=t_intensity;
    length=event.size();
    Create_dflux();
}

vector<double> Spectra::Filter_dflux(Histogram filter){

    double weight;
    vector<double> t_dflux;
    long dist = distance(filter.event.begin(),lower_bound(filter.event.begin(),filter.event.end(),event[0]));

    for(int i=0; i<length; i++){

        if(abs(filter.event[dist]-event[i]) > abs(filter.event[dist+1]-event[i])){dist+=1;}
        weight=filter.intensity[dist]+(filter.d_intensity[dist]) / (filter.d_event[dist]) * (filter.event[dist] - event[i]);

        if(filter.event[filter.length-1] - event[i]<0){weight=0;}
        if(filter.event[0] - event[i]>0){weight=0;}

        t_dflux.push_back(intensity[i]*(0.503)*event[i]*weight);

    }

    return t_dflux;
}

Spectra Spectra::Pass_filter(Histogram filter){

    dflux=Filter_dflux(filter);

    return *this;
}

void Spectra::Create_dflux(){
    for(int i=0; i<length; i++){
        dflux.push_back(intensity[i]*(0.503)*(event[i]));
        //0.503 is for assuming intensity is erg/s/cm^2/cm
        //and that event is in A
    }

    dflux.push_back(intensity[length-1]*(0.503)*event[length-1]);

    return;
}

void Spectra::Set_magnitude(double mag, double t_flux){

    magnitude=mag;

    double b=pow(10, -mag / 2.5)*v_zp / t_flux;

    for(int i=0; i<length+1; i++){
        dflux[i]=b*dflux[i];
        intensity[i]=b*intensity[i];
    }

    return;
}

double Spectra::Relative_flux(Histogram m_order, vector<Histogram> m_orders, double r_flux){

    double total_flux=0;
    double m_flux=Calc_flux(m_order); //Calc_flux(m_order)
    double a=0;

    for(int i=0; i<m_orders.size(); i++){

        a=Calc_flux(m_orders[i]);
        total_flux+=a;
    }

    return r_flux*(m_flux / total_flux);

}

double Spectra::Relative_flux(Histogram m_order, vector<Histogram> m_orders){

    return Calc_flux(m_order);

}

double Spectra::Calc_flux(const Histogram& filter){

   double flux = 0 ;
    vector<double> m_dflux= Filter_dflux(filter);


    for(int i=0; i<length; i++){

        flux+=(0.5)*(m_dflux[i+1]+m_dflux[i])*(d_event[i]);

    }

    return flux;

}

double Spectra::Calc_flux(){

    long double flux=0;

    for(int i=0; i<length; i++){

        flux=flux+(0.5)*(dflux[i+1]+dflux[i])*(d_event[i]);

    }

    magnitude=2.5*log10(v_zp / flux);
    return flux;

}

#endif //RANDOMGENERATOR_SPECTRA_H
