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
    std::vector<double> get_spectrum(std::vector<double> wavelength);

private:
    /*!
     * Integrates the \See Source::spectral_density() function between limits a and b.
     *
     * This is a simple integrator, which integrates the \see Source::spectral_density function between a and b.
     * It uses a simple aproximation by deviding the interval [a,b] in n parts and sum
     * \f[
     * I = \int_{a}^{b}(s(\lambda) d\lambda \approx \sum_{i=0}^{n} s(a + i*\frac{b-a}{n} * \frac{(b-a)}{n}
     * \f]
     * @param a lower wavelength limit
     * @param b upper wavelength limit
     * @param n number of subintervalls
     * @return integrated spectrum within [a, b]
     *
     * \todo This integrator should be replaved with a more accurate one. For highly unresolved spectra this integrator
     * might not be very precise.
     */
    double integral_s(double a, double b, int n);

};

/*!
 * Constant spectral density
 *
 * This class implements a constant spectral density
 */
class Constant : public Source{
public:
    /* Constructor */
    Constant();
    /*!
     * Constructor
     * @param value constant spectral density value
     */
    Constant(double value);
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
#endif // SOURCE_H
