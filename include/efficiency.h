#ifndef EFFICIENCY_H
#define EFFICIENCY_H

#include <string>
#include <vector>
#include <map>

class Efficiency {
public:
    Efficiency();

    virtual ~Efficiency();

    virtual std::vector<double> get_efficiency(int order, std::vector<double> &wavelength);

    virtual std::vector<double> get_efficiency(int order, std::vector<double> &wavelength, int N);

private:

};

class ConstantEfficiency : public Efficiency {
public:
    ConstantEfficiency(double efficiency);

    std::vector<double> get_efficiency(int order, std::vector<double> &wavelength);

    std::vector<double> get_efficiency(int order, std::vector<double> &wavelength, int N);

private:
    double eff;
};

/**
 * \class GratingEfficiency
 * \brief implements the efficiency curve of a echelle grating based on theory
 *
 */
class GratingEfficiency : public Efficiency {
public:
    /**
     * Constructor.
     * @param peak_efficiency peak efficiency of the grating
     * @param alpha alpha angle
     * @param blaze blaze angle
     * @param gpmm grooves per mm
     */
    GratingEfficiency(double peak_efficiency, double alpha, double blaze, double gpmm);

    std::vector<double> get_efficiency(int order, std::vector<double> &wavelength);

    std::vector<double> get_efficiency(int order, std::vector<double> &wavelength, int N);

private:
    double peak_efficiency;
    double alpha;
    double blaze;
    double gpmm;

    double calc_eff(double scalingfactor, int order, double alpha, double blaze, double wl, double n);

};


class CSVEfficiency : public Efficiency {
public:
    CSVEfficiency(std::string path);

    std::vector<double> get_efficiency(int order, std::vector<double> &wavelength);

    std::vector<double> get_efficiency(int order, std::vector<double> &wavelength, int N);

private:
    std::vector<double> wl;
    std::vector<double> ef;
    std::map<double, double> data;
};

#endif // EFFICIENCY_H
