#include <cmath>
#include <algorithm>
#include "source.h"
#include <CCfits/CCfits>
#include <helper.h>
#include <source.h>
/// mark fmt as header only
#define FMT_HEADER_ONLY

#include <fmt/format.h>


// integration routine
template<typename Method, typename F, typename Float>
double integrate(F f, Float a, Float b, int steps, Method m) {
    double s = 0;
    double h = (b - a) / steps;
    for (int i = 0; i < steps; ++i)
        s += m(f, a + h * i, h);
    return h * s;
}

// methods
class rectangular {
public:
    enum position_type {
        left, middle, right
    };

    rectangular(position_type pos) : position(pos) {}

    template<typename F, typename Float>
    double operator()(F f, Float x, Float h) const {
        switch (position) {
            case left:
                return f(x);
            case middle:
                return f(x + h / 2);
            case right:
                return f(x + h);
        }
    }

private:
    const position_type position;
};

Source::Source() : integration_steps(10), shift(1.), list_like(true) {
}

Source::~Source() {

}

std::vector<double> Source::get_interpolated_spectral_density(std::vector<double> wavelength) {
    if (!this->list_like) {
        std::vector<double> spectrum;
        std::vector<double> diff;

        diff.reserve(wavelength.size());
        for (int i = 0; i < wavelength.size() - 2; i++)
            diff.push_back(this->shift * (wavelength[i + 1] - wavelength[i]));
        diff.push_back(diff.back());
        diff.push_back(diff.back());

        auto f = [this](double wl) { return this->get_spectral_density(wl); };

        for (int i = 0; i < wavelength.size(); i++) {
            double a = this->shift * wavelength[i] - diff[i] / 2.;
            double b = this->shift * wavelength[i] + diff[i] / 2.;
//            double result2 = this->integral_s(a, b, this->integration_steps) / diff[i];
            double result = s_val * integrate(f, a, b, integration_steps, rectangular(rectangular::middle)) / diff[i];
            spectrum.push_back(result);
        }

        return spectrum;
    } else {
        std::vector<double> spectrum;

        spectrum.reserve(wavelength.size());
        for (double wl : wavelength) {
            spectrum.push_back(get_spectral_density(wl));
        }

        return spectrum;
    }
}

std::vector<double> Source::get_photon_flux(std::vector<double> wavelength) {
    if (!this->list_like) {
        // convert microwatts / micrometer to photons / s per wavelength intervall
        double ch_factor = 5.03E12;
        std::vector<double> flux;
        std::vector<double> diff;

        for (int i = 0; i < wavelength.size() - 2; i++) {
            diff.push_back(wavelength[i + 1] - wavelength[i]);
        }
        diff.push_back(diff.back());
        diff.push_back(diff.back());
        // expected to be microwatts/micrometer
        std::vector<double> spectrum = this->get_interpolated_spectral_density(wavelength);

        for (int i = 0; i < wavelength.size(); i++) {
            flux.push_back(spectrum[i] * wavelength[i] * ch_factor * diff[i]);
        }

        return flux;
    } else {
        std::vector<double> spectrum;

        for (int i = 0; i < wavelength.size(); i++) {
            // expected to be photons / s for list like sources
            spectrum.push_back(get_spectral_density(wavelength[i]));
        }

        return spectrum;
    }
}

void Source::set_integration_steps(int n) {
    this->integration_steps = n;
}

double Source::integral_s(double a, double b, int n) {
    double step = (b - a) / n;  // width of each small rectangle
    double area = 0.0;  // signed area
    for (int i = 0; i < n; i++) {
        area += this->get_spectral_density(a + (i + 0.5) * step) * step; // sum up each small rectangle
    }
    return area;
}

void Source::set_doppler_shift(double shift) {
    this->shift = 1. + shift / 2.99792458E8;
}


double Constant::get_spectral_density(double wavelength) {
    return this->value;
}

Constant::Constant(double value) : value(value) {
    list_like = false;
    name = fmt::format("Constant {0}", value);
}

Constant::Constant() {
    this->value = 0.01; // uW per um (micro watts per micro meter)
}

Blackbody::Blackbody(double T, double magnitude, double telescope_area) : StellarSource(magnitude, telescope_area),
                                                                          T(T) {
    this->list_like = false;
    name = fmt::format("Blackbody T: {0}, mag: {1}, telescope: {2}", T, magnitude, telescope_area);
    min_w = 0; // Maybe set a cutoff based off intensity?
    max_w = 10;

    calc_flux_scale();

}

double Blackbody::planck(const double &T, const double &wavelength) {
    const double hPlanck = 6.62606896e-34; //J * s;
    const double speedOfLight = 2.99792458e8; //m / s;
    const double kBoltzmann = 1.3806504e-23; //J / K
    double a = 2.0 * hPlanck * speedOfLight * speedOfLight;
    double b = hPlanck * speedOfLight / (wavelength * kBoltzmann * T);
    double c = 1E0; // Conversion factor for intensity
    double intensity = c * a / (pow(wavelength, 5) * (exp(b) - 1.0)); // (J / s ) / m^3 -> (uW) / ( m^2 * um )
    return intensity;
}

double Blackbody::get_spectral_density(double wavelength) {
    return this->planck(this->T, wavelength / 1E6);
}

PhoenixSpectrum::PhoenixSpectrum(std::string spectrum_file, std::string wavelength_file, double magnitude,
                                 double telescope_area) : StellarSource(magnitude, telescope_area) {
    list_like = false;
    name = fmt::format("Phoenix file: {0}, mag: {1}, telescope: {2}", spectrum_file, magnitude, telescope_area);
    this->read_spectrum(spectrum_file, wavelength_file);
    calc_flux_scale();
}

void PhoenixSpectrum::read_spectrum(std::string spectrum_file, std::string wavelength_file) {
// read in wavelength file
    std::unique_ptr<CCfits::FITS> ptr_FITS_wl(new CCfits::FITS(wavelength_file, CCfits::Read, true));
    CCfits::PHDU &wl = ptr_FITS_wl->pHDU();
    std::valarray<double> contents_wl;
    wl.readAllKeys();
    wl.read(contents_wl);

    std::unique_ptr<CCfits::FITS> ptr_FITS_spectrum(new CCfits::FITS(spectrum_file, CCfits::Read, true));
    CCfits::PHDU &spec = ptr_FITS_spectrum->pHDU();
    std::valarray<double> contents_spec;
    spec.readAllKeys();
    spec.read(contents_spec);
    // find wavelength limits
    // convert contents_spec from erg/s/cm^2/cm to uW/m^2/um
    for (long j = 0; j < contents_wl.size() - 1; ++j) {
        this->data[contents_wl[j] / 10000.] = 0.1 * contents_spec[j];
    }
}

double PhoenixSpectrum::get_spectral_density(double wavelength) {
    return interpolate(this->data, wavelength);
}

CoehloSpectrum::CoehloSpectrum(std::string spectrum_file, double magnitude, double telescope_area)
        : StellarSource(magnitude, telescope_area) {
    this->read_spectrum(std::move(spectrum_file));
    list_like = false;
    name = fmt::format("Coehlo spectrum: {0}, mag: {1}, telescope: {2}", spectrum_file, magnitude, telescope_area);
    calc_flux_scale();
}

void CoehloSpectrum::read_spectrum(std::string spectrum_file) {
//// read in wavelength file
    std::unique_ptr<CCfits::FITS> ptr_FITS_spectrum(new CCfits::FITS(spectrum_file, CCfits::Read, true));
    CCfits::PHDU &spec = ptr_FITS_spectrum->pHDU();
    std::valarray<double> contents_spec;
    spec.readAllKeys();
    spec.read(contents_spec);
    long ax1(spec.axis(0));
    long max_idx = ax1;
    long min_idx = 0;
    std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(spectrum_file, CCfits::Read));
    double minwl;
    pInfile->pHDU().readKey("CRVAL1", minwl);

    double wld;
    pInfile->pHDU().readKey("CDELT1", wld);

    // find wavelength limits
    // convert contents_spec from erg/s/cm^2/A to uW/s/m^2/um
    for (long j = min_idx; j < max_idx; ++j) {
        this->data[(minwl + j * wld) / 10000.] = (1E7) * contents_spec[j];
    }

}

double CoehloSpectrum::get_spectral_density(double wavelength) {

    return interpolate(this->data, wavelength);

}

CustomSpectrum::CustomSpectrum(double magnitude, double telescope_area, const std::string spectrum_file,
                               std::string wave_file) : StellarSource(magnitude, telescope_area) {
    list_like = false;
    name = fmt::format("custom : {0}, mag: {1}, telescope: {2}", spectrum_file, magnitude, telescope_area);
    read_spectrum(spectrum_file, wave_file);
    calc_flux_scale();

}

void CustomSpectrum::read_spectrum(const std::string spectrum_file, std::string wave_file) {
    //// read in wavelength file
    std::unique_ptr<CCfits::FITS> ptr_FITS_wl(new CCfits::FITS(wave_file, CCfits::Read, true));
    CCfits::PHDU &wl = ptr_FITS_wl->pHDU();
    std::valarray<double> contents_wl;
    wl.readAllKeys();
    wl.read(contents_wl);

    long ax1(wl.axis(0));
    long max_idx = ax1;
    long min_idx = 0;
    double min_wavelength_angstrom = min_w * 10000.;
    double max_wavelength_angstrom = max_w * 10000.;

    // find wavelength limits
    for (long j = 0; j < ax1; ++j) {
        if (contents_wl[j] < min_wavelength_angstrom)
            min_idx = j;
        if (contents_wl[j] > max_wavelength_angstrom)
            max_idx = j;
    }

    std::unique_ptr<CCfits::FITS> ptr_FITS_spectrum(new CCfits::FITS(spectrum_file, CCfits::Read, true));
    CCfits::PHDU &spec = ptr_FITS_spectrum->pHDU();
    std::valarray<double> contents_spec;
    spec.readAllKeys();
    spec.read(contents_spec);
    std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(spectrum_file, CCfits::Read));

    // find wavelength limits
    // convert contents_spec from erg/s/cm^2/A to uW/s/m^2/um
    for (long j = min_idx; j < max_idx; ++j) {
        this->data[contents_wl[j] / 10000.] = (1E7) * contents_spec[j];
    }
}

double CustomSpectrum::get_spectral_density(double wavelength) {
    return interpolate(this->data, wavelength);
}

LineList::LineList(const std::string linelist_file) {
    name = fmt::format("LineList file: {0}", linelist_file);
    list_like = true;
    this->read_spectrum(linelist_file);
}

void LineList::read_spectrum(std::string linelist_file) {
    std::ifstream file(linelist_file.c_str());
    for (CSVIterator loop(file); loop != CSVIterator(); ++loop) {
        std::cout << (*loop)[0] << std::endl;
        event.push_back(stod((*loop)[0]));
        intensity.push_back(stod((*loop)[1]));
        this->data.insert(std::pair<double, double>(stod((*loop)[0]), stod((*loop)[1])));
    }
}

std::vector<double> LineList::get_interpolated_spectral_density(std::vector<double> wavelength) {
    std::vector<double> spectrum;
    for (auto const &wl: wavelength) {
        spectrum.push_back(double(data[wl]));
    }
    return spectrum;
}

double LineList::get_spectral_density(double wavelength) {
    return data[wavelength];
}

std::vector<double> LineList::get_wavelength() {
    std::vector<double> wavelength;
    for (auto m: this->data) {
        wavelength.push_back(m.first * shift);
    }
    return wavelength;
}

CalibrationSource::CalibrationSource() {
    stellar_source = false;
}

StellarSource::StellarSource(double magnitude, double telescope_area) : mag(magnitude), telescope_area(telescope_area) {

}

void StellarSource::calc_flux_scale() {
    // V-band filter
    std::vector<double> v_filter_wl{0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6,
                                    0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7};
    std::vector<double> v_filter_tp{0, 0.03, 0.163, 0.458, 0.78, 0.967, 1, 0.973, 0.898, 0.792, 0.684, 0.574, 0.461,
                                    0.359, 0.27, 0.197, 0.135, 0.081, 0.045, 0.025, 0.017, 0.013, 0.009, 0};
    std::map<double, double> v_filter;
    for (size_t i = 0; i < v_filter_wl.size(); ++i)
        v_filter[v_filter_wl[i]] = v_filter_tp[i];
    s_val = 1.0;

    // get total flux in filter range
    std::vector<double> wl;
    double lower_wl_limit = std::max(min_w, v_filter_wl.front());
    double upper_wl_limit = std::min(max_w, v_filter_wl.back());

    int n = 1000000;
    double step = (upper_wl_limit - lower_wl_limit) / n;
    for (int i = 0; i < n; ++i)
        wl.push_back(lower_wl_limit + step * i);

    std::vector<double> spec = this->get_interpolated_spectral_density(wl);

    // get flux * V-filter
    double total_flux = 0.;
    for (int i = 0; i < n; i++) {
        total_flux += spec[i] * interpolate(v_filter, lower_wl_limit + step * i) * step;
    }

    s_val = pow(10, mag / (-2.5)) * v_zp / total_flux * telescope_area;
}


IdealEtalon::IdealEtalon(double d, double n, double theta, double R, double I) : d(d / 1000.), n(n), theta(theta),
                                                                                 R(R), I(I) {
    this->cF = this->coefficient_of_finesse(R);
    this->integration_steps = 10;
    this->list_like = false;
}

double IdealEtalon::coefficient_of_finesse(double R) {
    return 4. * R / ((1. - R) * (1. - R));
}

double IdealEtalon::T(double wl, double theta, double d, double n, double cF) {
    //delta = (2. * math.pi / wl) * 2. * n * math.cos(theta) * d
    //return 1. / (1. + cF * np.sin(0.5 * delta) ** 2)

    double delta = (2. * M_PI * n * cos(theta) * d) / wl;
    double sind = sin(delta);
    return 1. / (1. + cF * sind * sind);
}

double IdealEtalon::get_local_efficiency(double wavelength) {
    return this->T(wavelength / 1E6, theta, d, n, cF);
}

double IdealEtalon::get_spectral_density(double wavelength) {
    return I * get_local_efficiency(wavelength);
}
