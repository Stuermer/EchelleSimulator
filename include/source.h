#ifndef SOURCE_H
#define SOURCE_H

#include <vector>
#include <map>
/*!
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

    /*!
     * \overload virtual std::vector<double> get_spectral_density(std::vector<double> wavelength);
     * Returns the spectral density of the Source at the given wavelength vector.
     *
     * @param wavelength wavelength vector
     * @return spectral density at given wavelength
     */
    virtual std::vector<double> get_spectral_density(std::vector<double> wavelength);

    /*!
     * \fn virtual double get_spectral_density(double wavelength)
     * This function returns the spectral density at a given wavelength. It is the essential function for all subclasses.
     *
     * \see Source::get_spectrum() will use this function to integrate over it to retreive a spectrum for a given wavelength vector.
     * @param wavelength wavelength
     * @return spectral density
     */
    virtual double get_spectral_density(double wavelength);

    /*!
     * Returns spectrum at given wavelength
     *
     * This function returns the integrated spectral density for a given wavelength vector.
     * @param wavelength wavelength vector
     * @return spectrum at given wavelength
     */
    virtual std::vector<float> get_spectrum(std::vector<double> wavelength);

    /*!
     * Applies a spectral shift on the spectrum to simulate radial velocity shifts.
     *
     * @param shift doppler shift in [m/s]
     */
    void set_doppler_shift(double shift);

    /*!
     * Sets the number of sub steps of the integrator.
     * @param n number of subintervalls
     */
    void set_integration_steps(int n);

    /*!
     * Scales the spectral density of the source by converting to photon density and normalizing against integrated photon flux
     *
     * Spectral density is assumed to be in the units [micro watt] / ([micro meter] * [meter]^2). To convert to photon density we divide the spectral density by
     * the energy in a photon at a specific wavelength (Planck's Equation). This results in the multiplied factors Source::inten_pho and wavelength
     * which is represented by (a + (i + 0.5) * step). We integrate over all available wavelengths, Source::min_w to Source::max_w, to obtain
     * the photon flux. We then produce a scaling factor Source::s_val for the spectral density by comparing the photon flux
     * to the flux from the star Vega. We then normalize so that our source is at a fixed magnitude Source::mag with respect to Vega.
     */
    void scale_spectral_density();

protected:

    /// Source apparent magnitude
    double mag;
    /// Reference flux obtained from integration of vega over bessel filter (units are photons/m^2/s)
    double v_zp=8660006000.0;
    /// Scaling factor used in Source::scale_spectral_density() for normalization of source spectral_density against Source::v_zp
    double s_val  = 1.0;

    /// minimum wavelength recorded for source [micro meters]
    double min_w = 0.45;
    /// maximum wavelength recorded for source [micro meters]
    double max_w = 0.85;

private:

    /*!
     * Integrates the \see{Source::spectral_density()} function between limits a and b.
     *
     * This is a simple integrator, which integrates the \see{Source::spectral_density} function between a and b.
     * It uses a simple aproximation by deviding the interval [a,b] in n parts and sum
     * \f[
     * I = \int_{a}^{b}(s(\lambda) d\lambda \approx \sum_{i=0}^{n} s(a + (i+0.5)frac{b-a}{n}) * \frac{(b-a)}{n}
     * \f]
     * @param min_w lower wavelength limit
     * @param max_w upper wavelength limit
     * @param n number of subintervalls
     * @return integrated spectrum within [a, b]
     *
     * \todo This integrator should be replaved with a more accurate one. For highly unresolved spectra this integrator
     * might not be very precise.
     */
    double integral_s(double a, double b, int n);

    double shift; ///< current doppler shift
    int integration_steps; ///< number of steps for the integrator
};

/*!
 * \class Constant
 * \brief Implements constant spectral density.
 *
 * This class implements a constant spectral density.
 * \f[
 * s(\lambda) = const.
 * \f]
 */
class Constant : public Source{
public:
    /* Constructor */
    Constant();
    /*!
     * Constructor with a wavelength range of 0 to 1 meter
     * @param value constant spectral density value
     */
    Constant(double value, double min_w, double max_w);
    /*!
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

/*!
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
class IdealEtalon : public Source{
public:
    /*!
     * Constructor.
     *
     * @param d mirror distance in mm
     * @param n refractive index of the medium between mirrors
     * @param theta angle of incidence
     * @param R reflectivity of the mirrors
     */
    IdealEtalon(double d, double n, double theta, double R);
    /*!
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

    double FSR();

    double F();

    /*!
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
    /*!
     * Spectral density at given wavelegnth.
     *
     * @param wavelength wavelength [micron]
     * @return Spectral density at given wavelength
     */
    double get_spectral_density(double wavelength);
private:
    double d;
    double n;
    double theta;
    double R;
    double cF;
};


/*!
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
class Blackbody : public Source{
public:
    /*!
     * Constructor
     * @param T Temperature [K]
     */
    Blackbody(double T, double mag);
    /*!
     * Planck function for spectral density of a blackbody with Temperature T
     * @param T Temperature [K]
     * @param wavelength wavelength [m]
     * @return spectral density
     */
    double planck(const double& T, const double& wavelength);
    /*!
     * spectral density of a blackbody
     * \f[
     * s(\lambda) = \frac{2hc^2}{\lambda^5}\frac{1}{\exp{\frac{hc}{\lambda k_B T}}-1}
     * \f]
     * @param wavelength wavelength [micro meters]
     * @return spectral density of a blackbody at given wavelength [micro Watt] / ([micro meter] * [meter]^2 )
     */
    double get_spectral_density(double wavelength);

private:
    double T; ///< Temperature [K]
};

class PhoenixSpectrum : public Source{
public:
    PhoenixSpectrum(std::string spectrum_file, std::string wavelength_file, double mag);
    /*!
     * read in the spectrum from file between min_wavelength and max_wavelength
     * @param spectrum_file file path of spectrum
     * @param wavelength_file file path of wavelength
     * @param min_wavelength minimum wavelength in micrometer
     * @param max_wavelength maximum wavelength in micrometer
     */
    void read_spectrum(std::string spectrum_file, std::string wavelength_file, double mag);
    double get_spectral_density(double wavelength);

private:

    std::map<double, double> data;

};

/*!
 * \class LineList
 * \brief Implements line list spectrum
 *
 */
class LineList : public Source{
public:
    LineList(std::string linelist);
    void read_spectrum(std::string linelist);
    double get_spectral_density(double wavelength);
    std::vector<float> get_spectrum(std::vector<double> wavelength);
    std::vector<double> get_wavelength();
private:
    std::map<double, double> data;
};
#endif // SOURCE_H
