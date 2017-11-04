//
// Created by zrobertson on 6/27/17.
//

#ifndef RANDOMGENERATOR_SPECTRA_H
#define RANDOMGENERATOR_SPECTRA_H

#include "Histogram.h"
#include <iostream>

/*!
 * \class Spectra
 * \brief Does computations with sources (such as sampling)
 *
 * This class allows for a variety of computations to be done with sources. It's basic
 * functionality is to allow filters to be applied to sources and to allow conversion between
 * energy flux and photon flux.
 *
 */

class Spectra:public Histogram{

public:
    Spectra();
    Spectra(string path);
    Spectra(string path_1, string path_2);

    Spectra(vector<double> m_event, vector<double> m_intensity);
    Spectra(Spectra spectra, Histogram filter);
    //Filters the spectra and adds the result as an object
    //Technically deprecated because of Pass_filter

    /// The reference flux units :: photons/m^2/s (obtained by integrating vega over a bessel filter)
    double v_zp=8660006000.0; //The reference flux is obtained by integrating vega
                                // over a bessel filter and has units photons/m^2/s

    double magnitude;

    /// Photon flux units:: photons/s/cm^2/A
    vector<double> dflux; //Units are flux per wavelength  photons/s/cm^2/A

    /*!
     * Creates a photon flux vector from the intesnity vector. Uses the formula
     * Flux_p = Intensity * ch_factor * lambda where ch_factor = 5.03*10^10 uW/m^2/um
     * with the assumption that lambda is the wavelength in units :: um
     *
     */
    void Create_dflux();
    void Set_magnitude(double mag, double t_flux);

    /// Passes a filter over the Spectra and recomputes @dflux
    Spectra Pass_filter(Histogram filter); //Passes the filter over Spectra
    /// The function that recomputes @dflux according to a filter
    vector<double> Filter_dflux(Histogram filter);

    /*!
     * Returns the flux of single spectra assuming that it is one of several
     * spectra that have a normalized total flux.
     *
     * @param m_order the reference histogram
     * @param m_orders the list of histograms to be referenced
     * @param r_flux the assumed total flux of all the orders
     * @return
     */
    double Relative_flux(Histogram m_order, vector<Histogram> m_orders, double r_flux);
    //Calculates flux through m_order assuming
    // there are m_orders and that r_flux is conserved (slower by a factor of num(m_orders))

    /// Very similar to the other method but relaxes the normalization assumption
    double Relative_flux(Histogram m_order, vector<Histogram> m_orders); //Same as above but assumes r_flux = total_flux

    double Calc_flux(const Histogram& filter);
    double Calc_flux();

    /*!
     * Simulates the number of photons that a telescope would collect from a spectra
     * over a certain time. (Not currently implemented)
     *
     * @param m_time units :: seconds
     * @param m_area units :: m^2
     * @return
     */
    double Photon_count(double m_time, double m_area);
    double Exp_time(double s_noise);

};

Spectra::Spectra():Histogram(){

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
    //0.503 is for assuming intensity is erg/s/cm^2/cm
    //and that event~(wavelength) is in A

    //5.03*10^10 is for assuming intensity is uW/m^2/um
    //and that event is in um

    double ch_factor = 5.03E10;

    for(int i=0; i<length; i++){
        dflux.push_back(intensity[i]*(ch_factor)*(event[i]));

    }

    dflux.push_back(intensity[length-1]*(ch_factor)*event[length-1]);

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

        long double flux = 0;

        for (int i = 0; i < length; i++) {

            flux = flux + (0.5) * (dflux[i + 1] + dflux[i]) * (d_event[i]);

        }

        magnitude = 2.5 * log10(v_zp / flux);
        return flux;

}

#endif //RANDOMGENERATOR_SPECTRA_H
