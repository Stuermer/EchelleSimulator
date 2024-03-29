#include <vector>
#include <iostream>
#include "PSF.h"
#include "H5Cpp.h"
#include "helper.h"


herr_t
dataset_info(hid_t loc_id, const char *name, const H5L_info_t *linfo, void *opdata) {
    hid_t ds;
    auto dataset_names = reinterpret_cast< std::vector<std::string> * >(opdata);
    ds = H5Dopen2(loc_id, name, H5P_DEFAULT);

    dataset_names->push_back(name);
    // std::cout << "Name : " << name << std::endl;
    H5Dclose(ds);
    return 0;
}

PSF::PSF() = default;

PSF::~PSF() = default;

PSF_ZEMAX::PSF_ZEMAX(const std::string& filename, int fiber_number) {

    auto *file = new H5::H5File(filename, H5F_ACC_RDONLY);
    std::cout << std::endl << "Iterating over elements in the file" << std::endl;
    auto *rootGr = new H5::Group(file->openGroup("fiber_" + std::to_string(fiber_number)));
    std::vector<std::string> group_names;
    // H5::Group *rootGr = new H5::Group (file->openGroup("/"));
    herr_t idx = H5Literate(rootGr->getId(), H5_INDEX_NAME, H5_ITER_INC, NULL, file_info, &group_names);

    for (auto &gn : group_names) {
        if (gn.find("psf") != std::string::npos) {
            // group names are psf_order_100    get order number here
            int order = std::stoi(gn.substr(10));
            std::vector<PSFdata> psf_vec;
            std::vector<double> pixelsampling;

            std::vector<std::string> wl_names;
            const std::string& dt_path = gn;
            auto *wlGr = new H5::Group(rootGr->openGroup(dt_path));


            herr_t idx2 = H5Literate(wlGr->getId(), H5_INDEX_NAME, H5_ITER_INC, nullptr, dataset_info, &wl_names);
            for (auto &w : wl_names) {
                H5::DataSet dataset = rootGr->openDataSet(dt_path + "/" + w);

                auto *attr = new H5::Attribute(dataset.openAttribute("wavelength"));
                auto *type = new H5::DataType(attr->getDataType());
                double wavelength = 0.;
                attr->read(*type, &wavelength);

                *attr = H5::Attribute(dataset.openAttribute("dataSpacing"));
                *type = H5::DataType(attr->getDataType());
                double read_sampling = 0.;
                attr->read(*type, &read_sampling);
                this->pixelsampling = read_sampling;


                H5::DataSpace dspace = dataset.getSpace();
                int ndims = dspace.getSimpleExtentNdims();
                hsize_t dims[ndims];

                dspace.getSimpleExtentDims(dims, NULL);

                auto *buffer = new double[(unsigned long) dims[0] * (unsigned long) dims[1]];

                H5::DataSpace mspace1 = H5::DataSpace(2, dims);
                dataset.read(buffer, H5::PredType::NATIVE_DOUBLE, mspace1, dspace);


                Matrix raw_psf((unsigned long) dims[0], (unsigned long) dims[1]);

                for (int i = 0; i < (unsigned long)  dims[0]; ++i) {
                    for (int j = 0; j < (unsigned long)  dims[1]; ++j) {
                        raw_psf.data[i][j] = (float) buffer[i * (unsigned long)  dims[1] + j];
                    }
                }
                PSFdata p(wavelength, raw_psf);
                if (this->psfs.find(order) == this->psfs.end()) {
                    // not found
                    psf_vec.emplace_back(p);
                    this->psfs.insert(std::make_pair(order, psf_vec));
                } else {
                    this->psfs[order].emplace_back(p);
                }
            }
            delete wlGr;

            std::sort(this->psfs[order].begin(), this->psfs[order].end());
        }
    }
    delete rootGr;
    delete file;
}

Matrix PSF_ZEMAX::get_PSF(int order, double wavelength) {
    std::vector<double> distance;
    for (auto &psf : this->psfs[order]) {
        distance.push_back(psf.wavelength - wavelength);
    }
    std::vector<size_t> idx;
    idx = compute_sort_order(distance);
    return this->interpolate_PSF(this->psfs[order][idx[0]].psf, this->psfs[order][idx[1]].psf,
                                 this->psfs[order][idx[0]].wavelength, this->psfs[order][idx[1]].wavelength,
                                 wavelength);
}

Matrix PSF_ZEMAX::get_PSF_nocut(int order, double wavelength) {
    std::vector<double> distance;
    for (auto &psf : this->psfs[order]) {
        distance.push_back(psf.wavelength - wavelength);
    }
    std::vector<size_t> idx;
    idx = compute_sort_order(distance);
    return this->interpolate_PSF_nocut(this->psfs[order][idx[0]].psf, this->psfs[order][idx[1]].psf,
                                       this->psfs[order][idx[0]].wavelength, this->psfs[order][idx[1]].wavelength,
                                       wavelength);
}

Matrix PSF_ZEMAX::interpolate_PSF(Matrix * psf1, Matrix * psf2, double w1, double w2, double w) {
    double p1 = fabs((w - w1) / (w2 - w1));
    double p2 = fabs((w - w2) / (w2 - w1));
    double p_sum = p1 + p2;
    p1 /= p_sum;
    p2 /= p_sum;
    Matrix comb_psf(psf1->rows, psf1->cols);
    double psf1_total = psf1->sum();
    double psf2_total = psf2->sum();

    int minX, maxX, minY, maxY;
    minX = comb_psf.cols;
    minY = comb_psf.rows;
    maxX = 0;
    maxY = 0;

    // interpolate and cut PSF to smallest rectangle that contains values > 0.001.
    // When random values are drawn it helps to quicker find valid XY positions since we use the rejection method there
    float max_val = 0;
    for (int i = 0; i < comb_psf.rows; ++i) {
        for (int j = 0; j < comb_psf.cols; ++j) {
            double val = (p2 * psf1->data[i][j] / psf1_total + p1 * psf2->data[i][j] / psf2_total);
            if (val > 0.001) {
                minX = j < minX ? j : minX;
                minY = i < minY ? i : minY;
                maxX = j > maxX ? j : maxX;
                maxY = i > maxY ? i : maxY;
            }
            comb_psf.data[i][j] = val;
            if (val > max_val)
                max_val = val;
        }
    }

    int cenX = psf1->cols / 2;
    int cenY = psf1->rows / 2;


    int size_x = abs(((cenX - minX) > (maxX - cenX)) ? cenX - minX : maxX - cenX);
    int size_y = abs(((cenY - minY) > (maxY - cenY)) ? cenY - minY : maxY - cenY);
    Matrix result(size_y * 2, size_x * 2);

//    comb_psf.delete_n_rows_symmetrically(comb_psf.rows-size_y*2);
//    comb_psf.delete_n_cols_symmetrically(comb_psf.cols-size_x*2);
    //normalize maximum value to 1:
    for (int i = 0; i < result.rows; ++i) {
        for (int j = 0; j < result.cols; ++j) {
            result.data[i][j] = comb_psf.data[comb_psf.rows - size_y * 2 + i][comb_psf.cols - size_x * 2 + j] / max_val;
        }
    }
    return result;

}

Matrix PSF_ZEMAX::interpolate_PSF_nocut(Matrix * psf1, Matrix * psf2, double w1, double w2, double w) {
    double p1 = fabs((w - w1) / (w2 - w1));
    double p2 = fabs((w - w2) / (w2 - w1));
    double p_sum = p1 + p2;
    p1 /= p_sum;
    p2 /= p_sum;
    Matrix comb_psf(psf1->rows, psf1->cols);
    double psf1_total = psf1->sum();
    double psf2_total = psf2->sum();


    for (int i = 0; i < comb_psf.rows; ++i) {
//        const double * ptr1 = psf1.ptr<double>(i);
//        const double * ptr2 = psf2.ptr<double>(i);
//        double * ptr3 = comb_psf.ptr<double>(i);
        for (int j = 0; j < comb_psf.cols; ++j) {
            float val = (p2 * psf1->data[i][j] / psf1_total + p1 * psf2->data[i][j] / psf2_total);
            comb_psf.data[i][j] = val;
        }
    }

    return comb_psf;
}

PSF_gaussian::PSF_gaussian(double sigma, double aperture) : sigma(sigma) {
    this->ksize = ceil(sigma * aperture) * 2 + 1;
}


Matrix PSF_gaussian::get_PSF(int order, double wavelength) {

    return Matrix(30, 30);
}