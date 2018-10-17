#ifndef SOURCE_H
#define SOURCE_H

#include <vector>
#include <map>
/**
 * \class Source
 * \brief Base class of all spectral sources.
 *
 * This class is the base class of all spectral sources. Its purpose is to provide a common interface for all sources.
 * For implementing a new spectral source, inherit from this class and implement the Source::get_spectral_density function.
 */
class Source
{
public:
    /* Constructor */
    Source();

    /* Destructor - it's important that this is virtual ! */
    virtual ~Source();

    /**
     * Returns spectrum at given wavelength
     *
     * This function returns the integrated spectral density for a given wavelength vector.
     * @param wavelength wavelength vector
     * @return spectrum at given wavelength
     */
    virtual std::vector<double> get_interpolated_spectral_density(std::vector<double> wavelength);

    /**
     * Applies a spectral shift on the spectrum to simulate radial velocity shifts.
     *
     * @param shift doppler shift in [m/s]
     */
    void set_doppler_shift(double shift);

    /**
     * Sets the number of sub steps of the integrator.
     * @param n number of sub-intervals
     */
    void set_integration_steps(int n);

    virtual std::vector<double> get_photon_flux(std::vector<double> wavelength);

    virtual std::vector<double> get_wavelength() {};

    bool is_list_like() { return list_like;};
    bool is_stellar_source() {return stellar_source; };
    std::string get_source_name() {return name;};

protected:
        /**
     * \fn virtual double get_spectral_density(double wavelength)
     * This function returns the spectral density at a given wavelength. It is the essential function for all subclasses.
     *
     * \see Source::get_spectrum() will use this function to integrate over it to retrieve a spectrum for a given wavelength vector.
     * @param wavelength wavelength
     * @return spectral density
     */
    virtual double get_spectral_density(double wavelength) {};

    /// whether the source is list-like or not
    bool list_like;
    /// stellar source or not?
    bool stellar_source;
    /// name of the source
    std::string name;
    /// current doppler shift
    double shift;
    /// Scaling factor used in Source::scale_spectral_density() for normalization of source spectral_density against Source::v_zp
    double s_val  = 1.0;

    /**
     * Integrates the \see{Source::spectral_density()} function between limits a and b.
     *
     * This is a simple integrator, which integrates the \see{Source::spectral_density} function between a and b.
     * It uses a simple aproximation by deviding the interval [a,b] in n parts and sum
     * \f[
     * I = \int_{a}^{b}(s(\lambda) d\lambda \approx \sum_{i=0}^{n} s(a + (i+0.5)frac{b-a}{n}) * \frac{(b-a)}{n}
     * \f]
     * @param n number of sub-intervals
     * @param a lower limit for integration [micron]
     * @param b upper limit for integration [micron]
     * @return integrated spectrum within [a, b]
     *
     * \todo This integrator should be replaced with a more accurate one. For highly unresolved spectra this integrator
     * might not be very precise.
     */
    double integral_s(double a, double b, int n);

    ///< number of steps for the integrator
    int integration_steps;
};

class CalibrationSource : public Source{
public:
    CalibrationSource();
};

class StellarSource : public Source{
public:
    StellarSource(double magnitude, double telescope_area);

    // std::vector<double> get_interpolated_spectral_density(std::vector<double> wavelength);

protected:
    /// Source apparent magnitude
    double mag;
    /// telescope collecting light area
    double telescope_area;
    /// Reference flux obtained from integration of vega over bessel filter (units are microwatts/m^2*micrometer)
    double v_zp = 3.68E-02;

    /// minimum wavelength recorded for source [micro meters]
    double min_w = 0;
    /// maximum wavelength recorded for source [micro meters]
    double max_w = 10.;

    /**
    *
    * Scales the spectral density of the source by converting to photon density and normalizing against integrated photon flux

    * Spectral density is assumed to be in the units [micro watt] / ([micro meter] * [meter]^2). To convert to photon
    * density we divide the spectral density by the energy in a photon at a specific wavelength (Planck's Equation).
    * This results in the multiplied factors Source::inten_pho and wavelength
    * which is represented by (a + (i + 0.5) * step). We integrate over all available wavelengths,
    * Source::min_w to Source::max_w, to obtain the photon flux. We then produce a scaling factor
    * Source::s_val for the spectral density by comparing the photon flux to the flux from the star Vega.
    * We then normalize so that our source is at a fixed magnitude Source::mag with respect to Vega.
    */
    void calc_flux_scale();
};

/**
 * \class Constant
 * \brief Implements constant spectral density.
 *
 * This class implements a constant spectral density.
 * \f[
 * s(\lambda) = const.
 * \f]
 */
class Constant : public CalibrationSource{
public:
    /* Default Constructor */
    Constant();
    /**
     * Constructor with a wavelength range of 0 to 1 meter
     * @param value constant spectral density value
     */
    Constant(double value);
    /**
     * Returns constant spectral density.
     *
     * @param wavelength wavelength
     * @return constant spectral density value
     */
    double get_spectral_density(double wavelength);
private:
    /* Constant spectral density value */
    double value;
};

/**
 * \class IdealEtalon
 * \brief Implements the spectral density of an ideal fabry-perot etalon.
 *
 * An ideal Fabry-Perot etalon has a transmission function that only depends on the distance of the mirrors,
 * the angle of incidence, the reflectivity of the mirrors and the refractive index of the medium between the mirrors.
 * \f[
 * s(\lambda) = \frac{1}{cF sin(\frac{\delta} {2})^2}
 * \f]
 * (\see IdealEtalon::T())
 *
 * It produces a comb-like spectrum that is has equidistant peaks in frequency.
 *
 * \todo implement FSR(), F() and other static functions.
 *
 */
class IdealEtalon : public CalibrationSource{
public:
    /**
     * Constructor.
     *
     * @param d mirror distance in mm
     * @param n refractive index of the medium between mirrors
     * @param theta angle of incidence
     * @param R reflectivity of the mirrors
     * @param I flux density [microWatts]/[micrometer] of the underlying light source
     */
    IdealEtalon(double d, double n, double theta, double R, double I);
    /**
     * Calculates the coefficient of Finesse.
     *
     * \f[ cF = \frac{4 R}{(1-R)^2}
     * \f]
     *
     * \warning
     * This is not what is typically called the finesse of an etalon, but the coefficient of finesse.
     *
     * @param R mirror reflectivity
     * @return coefficienct of finesse
     */
    static double coefficient_of_finesse(double R);

    /**
     * Transmission function of an ideal etalon.
     * \f[ T(\lambda) = \frac{1}{cF sin(\frac{\delta} {2})^2}
     * \f], where
     * \f[ \delta \left( \frac{2\pi}{\lambda} \right) 2nlcos(\Theta)
     * \f],
     * is the phase difference and \f$ F \f$ is the coefficient of finesse:
     * \f[ cF = \frac{4R}{(1-R)^2}
     * \f]
     *
     * @param wl wavelength [micron]
     * @param theta angle of incidence [rad]
     * @param d mirror distance [mm]
     * @param n refractive index
     * @param cF coefficient of finesse \see IdealEtalon::coefficient_of_finesse
     * @return transmission at given wavelength
     */
    static double T(double wl, double theta, double d, double n, double cF);

    /**
     * Spectral density at given wavelegnth.
     *
     * @param wavelength wavelength [micron]
     * @return Spectral density at given wavelength
     */

protected:
    double get_spectral_density(double wavelength);

private:
    double get_local_efficiency(double wavelength);
    double d;
    double n;
    double theta;
    double R;
    double cF;
    double I;
};

/**
 * \class Blackbody
 * \brief Implements a *blackbody spectrum.*
 *
 * This class implements the spectrum of a blackbody of a certain Temperature.
 * \see{Blackbody::get_spectral_density()}
 *
 *
 *  \f[
 *  \fbox{
 *  \begin{tikzpicture}
 *  \begin{axis}[domain=300E-9:4000E-9, samples=100]
 *  \addplot[color=red]{2.0 *6.62606896E-34* 2.9979E8 * 2.9979E8 / (x^5)*1./(exp(6.62606896E-34*2.997E8/(x*1.3806504E-23*3500))-1)};
 *  \end{axis}
 *  \end{tikzpicture}
 *  }
 *  \f]
 *
 */
class Blackbody : public StellarSource{
public:
    /**
     * Constructor
     * @param T Temperature [K]
     * @param magnitude visual magnitude
     * @param telescope_area telescope light collecting area [m^2]
     */
    Blackbody(double T, double magnitude, double telescope_area);
    /**
     * Planck function for spectral density of a blackbody with Temperature T
     * @param T Temperature [K]
     * @param wavelength wavelength [m]
     * @return spectral density
     */
    double planck(const double& T, const double& wavelength);
    /**
     * spectral density of a blackbody
     * \f[
     * s(\lambda) = \frac{2hc^2}{\lambda^5}\frac{1}{\exp{\frac{hc}{\lambda k_B T}}-1}
     * \f]
     * @param wavelength wavelength [micro meters]
     * @return spectral density of a blackbody at given wavelength [micro Watt] / ([micro meter] * [meter]^2 )
     */
    double get_spectral_density(double wavelength);

private:
    /// Temperature [K]
    double T;
};

/**
 * \class PhoenixSpectrum
 * \brief Implements using data from the Phoenix data repository. \see also helper.csv::download_phoenix()
 *
 */
class PhoenixSpectrum : public StellarSource{
public:
    /**
     * Constructor. Reads in spectrum and wavelength files and scales it for given visual magnitude.
     * @param spectrum_file file path for spectrum
     * @param wavelength_file file path for wavelength file
     * @param magnitude magnitude of star in the V-band
     * @param telescope_area telescope light collecting area [m^2]
     */
    PhoenixSpectrum(std::string spectrum_file, std::string wavelength_file, double magnitude,
                    double telescope_area);

    double get_spectral_density(double wavelength);
private:
    /**
     * read in the spectrum from file between min_wavelength and max_wavelength
     * @param spectrum_file file path of spectrum
     * @param wavelength_file file path of wavelength
     */
    void read_spectrum(std::string spectrum_file, std::string wavelength_file);
    // contains wavelength, spectrum data
    std::map<double, double> data;
};

class CoehloSpectrum : public StellarSource{
public:
    CoehloSpectrum(std::string spectrum_file, double magnitude, double telescope_area);
    void read_spectrum(std::string spectrum_file);
    double get_spectral_density(double wavelength);

private:

    std::map<double, double> data;

};

class CustomSpectrum : public StellarSource{
public:
    CustomSpectrum(double magnitude, double telescope_area, const std::string spectrum_file,
                       std::string wave_file);

    void read_spectrum(const std::string spectrum_file, std::string wave_file);

    double get_spectral_density(double wavelength);

private:

    std::map<double, double> data;

};

/**
 * \class LineList
 * \brief Implements line list spectrum
 *
 */
class LineList : public CalibrationSource{
public:
    LineList(const std::string linelist_file);
    void read_spectrum(std::string linelist_file);
    double get_spectral_density(double wavelength);
    std::vector<double> get_interpolated_spectral_density(std::vector<double> wavelength);
    std::vector<double> get_wavelength();

    std::vector<double> event;
    std::vector<double> intensity;


private:
    std::map<double, double> data;
};
#endif // SOURCE_H
