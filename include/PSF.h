#ifndef ECHELLESIMULATOR_PSF_H
#define ECHELLESIMULATOR_PSF_H

#include <string>
#include <map>
#include <vector>

class Matrix {
public:
    Matrix(size_t rows, size_t cols) {
        this->cols = cols;
        this->rows = rows;
        data = std::vector<std::vector<float>>(rows, std::vector<float>(cols, 0.));
    }
    Matrix(const Matrix& m){
        this->cols = m.cols;
        this->rows = m.rows;
        this->data = m.data;
    }

    Matrix(std::vector<std::vector<float>> mat) {
        cols = mat[0].size();
        rows = mat.size();
        data = mat;
    }

    Matrix() {

    }

    ~Matrix() = default;

    std::vector<std::vector<float>> data;
    size_t cols{};
    size_t rows{};

    float sum() {
        float total = 0;
        for (int i = 0; i < data.size(); ++i) {
            for (int j = 0; j < data[i].size(); ++j) {
                total += data[i][j];
            }
        }
        return total;
    }

    void delete_n_rows_symmetrically(int n) {
        rows -= n;
        data.erase(data.begin(), data.begin() + n);
        data.erase(data.end() - n, data.end());
    }

    void delete_n_cols_symmetrically(int n) {
        cols -= n;
        for (unsigned i = 0; i < rows; ++i) {
            data[i].erase(data[i].begin(), data[i].begin() + n - 1);
            data[i].erase(data[i].end() - n + 1, data[i].end());
        }
    }

    Matrix &operator=(const Matrix &M) {
        this->cols = M.cols;
        this->rows = M.rows;
        this->data = M.data;
        return *this;
    }

};

struct PSFdata {
    double wavelength;
    Matrix * psf;

    PSFdata(double w, const Matrix& p) : wavelength(w) {
        psf = new Matrix(p);
    };

    bool operator<(const PSFdata &str) const {
        return (wavelength < str.wavelength);
    }
};

/*!
 * \class PSF
 * \brief Class handles point spread functions
 *
 * This class handles point spread functions (PSFs). It's basic functionality is to deliver a PSF as a 2d matrix for a
 * given wavelength and order of an echelle spectrograph.
 *
 * Typically the PSF of an echelle spectrograph will vary across the CCD depending on the wavelength, the echelle order,
 * and the illumination of the optics. The PSF might not be stable from target to target, as illumination might vary due
 * to different coupling conditions and imperfect scrambling of the fibers.
 *
 * To implement own PSF classes, inherit from PSF and overwrite get_PSF function
 */
class PSF {
public:
    PSF();

    virtual ~PSF();

    virtual Matrix get_PSF(int order, double wavelength) = 0;

    virtual Matrix get_PSF_nocut(int order, double wavelength) = 0;

    double pixelsampling{};
};

/**
 * \class PSF_ZEMAX
 * \brief class representing a PSF as computed by ZEMAX
 *
 * Zemax returns a 2D array with a calculated Huygens PSF. It is usually saved in the spectrograph model.
 */
class PSF_ZEMAX : public PSF {
public:
    /**
     * Constructor
     * @param filename filename of the spectrograph model
     * @param fiber_number fiber number for the PSF model
     */
    PSF_ZEMAX(const std::string& filename, int fiber_number);

    Matrix get_PSF(int order, double wavelength) override;

    Matrix get_PSF_nocut(int order, double wavelength) override;

private:
    /**
     * Interpolates between neighboring PSFs as a simple linear sum
     * @param psf1 first PSF
     * @param psf2 second PSF
     * @param w1 first PSF wavelength
     * @param w2 second PSF wavelength
     * @param w wavelength to use for calculating PSF
     * @return
     */
    Matrix interpolate_PSF(Matrix * psf1, Matrix * psf2, double w1, double w2, double w);

    Matrix interpolate_PSF_nocut(Matrix * psf1, Matrix * psf2, double w1, double w2, double w);

    std::map<int, std::vector<PSFdata> > psfs;

};

/**
 * \class PSF_gaussian
 * \brief Class representing a gaussian PSF
 *
 * This class handels gaussian like PSF functions.
 */
class PSF_gaussian : public PSF {
public:
    /**
     * Creator
     * @param sigma Sigma of the gaussian in (oversampled) pixel
     * @param aperture size of the gaussian kernel (in number of oversampled pixels)
     */
    PSF_gaussian(double sigma, double aperture = 3.);

    /**
     * Returns gaussian PSF, independently of order and wavelength
     * @param order echelle order
     * @param wavelength  wavelength in microns
     * @return
     */
    Matrix get_PSF(int order, double wavelength);

private:
    // sigma of gaussian in px
    double sigma;
    int ksize;

};

#endif //ECHELLESIMULATOR_PSF_H
