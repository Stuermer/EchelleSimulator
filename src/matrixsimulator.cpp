#include "matrixsimulator.h"
#include <iostream>
#include <fstream>
#include <iterator>
#include "Histogram.h"
#include "H5Cpp.h"
#include <cmath>
#include "helper.h"
#include <hdf5_hl.h>
#include <CCfits/CCfits>
#include <CCfits/FITS.h>
#include <string>
#include <fstream>
#include <chrono>
#include <utility>
#include <omp.h>
#include "random_generator.h"

#pragma omp declare reduction(vec_int_plus : std::vector<int> : \
                              std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<int>())) \
                    initializer(omp_priv(omp_orig.size(), 0))

typedef struct transformation_hdf {
    float rotation;
    float scale_x;
    float scale_y;
    float shear;
    float translation_x;
    float translation_y;
    float wavelength;
} transformation_hdf;

MatrixSimulator::MatrixSimulator(const std::string path, int fiber_number, bool keep_ccd) {
    this->load_spectrograph_model(path, fiber_number, keep_ccd);
}

void MatrixSimulator::load_spectrograph_model(const std::string path, int fiber_number, bool keep_ccd) {
    if (path.find(".hdf") != std::string::npos) {
        std::cout << "Spectrograph file found: " << path << std::endl;
    } else {
        std::cout << "Spectrograph file not found: " << path << std::endl;
    }
    const H5std_string &filename(path);
    this->efficiencies.clear();
    this->orders.clear();
    this->raw_transformations.clear();
    this->sim_wavelength.clear();
    this->sim_spectra.clear();
    this->sim_efficiencies.clear();
    this->sim_matrices.clear();

    // open file readonly
    auto *h5file = new H5::H5File(filename, H5F_ACC_RDONLY);

    // read in spectrograph information
    H5::Group *spec = new H5::Group(h5file->openGroup("Spectrograph"));
    H5::Attribute *attr = new H5::Attribute(spec->openAttribute("blaze"));
    H5::DataType *type = new H5::DataType(attr->getDataType());
    double blaze, gpmm = 0;
    attr->read(*type, &blaze);
    spec->openAttribute("gpmm").read(*type, &gpmm);
    this->spec_info.blaze = blaze;
    this->spec_info.gpmm = gpmm;

    // read in field information and transformations and PSFs
    H5::Group *fiber = new H5::Group(h5file->openGroup("fiber_" + std::to_string(fiber_number)));
    this->fiber_number = fiber_number;
    attr = new H5::Attribute(fiber->openAttribute("field_height"));
    type = new H5::DataType(attr->getDataType());
    double field_height, field_width = 0;
    int sampling_input_x;
    int number_of_points;

    attr->read(*type, &field_height);
    fiber->openAttribute("field_with").read(*type, &field_width);
    attr = new H5::Attribute(fiber->openAttribute("sampling_input_x"));
    type = new H5::DataType(attr->getDataType());
    fiber->openAttribute("sampling_input_x").read(*type, &sampling_input_x);
    fiber->openAttribute("MatricesPerOrder").read(*type, &number_of_points);

    this->slit = new Slit(field_width, field_height, sampling_input_x);

    std::vector<std::string> group_names;
    herr_t idx = H5Literate(fiber->getId(), H5_INDEX_NAME, H5_ITER_INC, nullptr, file_info, &group_names);
    size_t dst_size = sizeof(transformation_hdf);
    size_t dst_offset[7] = {HOFFSET(transformation_hdf, rotation),
                            HOFFSET(transformation_hdf, scale_x),
                            HOFFSET(transformation_hdf, scale_y),
                            HOFFSET(transformation_hdf, shear),
                            HOFFSET(transformation_hdf, translation_x),
                            HOFFSET(transformation_hdf, translation_y),
                            HOFFSET(transformation_hdf, wavelength)};
    size_t dst_sizes[7] = {sizeof(float), sizeof(float), sizeof(float), sizeof(float), sizeof(float), sizeof(float),
                           sizeof(float)};

    bool file_contains_psfs = false;

    for (const auto &gn: group_names) {
        if (gn.find("psf") == std::string::npos) {
            transformation_hdf data[number_of_points];
            std::string table_name = "fiber_" + std::to_string(fiber_number) + "/" + gn;
            int order = std::stoi(gn.substr(5));
            this->orders.push_back(order);
            herr_t read = H5TBread_table(h5file->getId(), table_name.c_str(), dst_size, dst_offset, dst_sizes, &data);
            for (int i = 0; i < number_of_points; ++i) {
                raw_transformation t;
                t.wavelength = data[i].wavelength;
                if (t.wavelength < this->wavelength_limit_min)
                    wavelength_limit_min = t.wavelength;
                if (t.wavelength > this->wavelength_limit_max)
                    wavelength_limit_max = t.wavelength;
                std::vector<float> result;
                result.push_back(data[i].scale_x);
                result.push_back(data[i].scale_y);
                result.push_back(wrap_rads(data[i].shear));
                result.push_back(data[i].rotation);
                result.push_back(data[i].translation_x);
                result.push_back(data[i].translation_y);
                t.decomposed_matrix = result;
                t.transformation_matrix = compose_matrix(result);
                t.order = order;
                this->raw_transformations[order].push_back(t);
            }
            std::sort(this->raw_transformations[order].begin(), this->raw_transformations[order].end(),
                      [](const raw_transformation &x, const raw_transformation &y) {
                          return (x.wavelength < y.wavelength);
                      });


        } else {
            file_contains_psfs = true;
        }
    }

    std::sort(this->orders.begin(), this->orders.end());
    this->min_order = *std::min_element(std::begin(this->orders), std::end(this->orders));
    this->max_order = *std::max_element(std::begin(this->orders), std::end(this->orders));
    this->n_orders = this->max_order - this->min_order;
    std::cout << "Read in " << this->n_orders << " Orders" << std::endl;
    this->calc_splines();

    if (!keep_ccd) {
        // read in CCD information
        H5::Group *ccd = new H5::Group(h5file->openGroup("CCD"));
        attr = new H5::Attribute(ccd->openAttribute("Nx"));
        type = new H5::DataType(attr->getDataType());
        int Nx = 0;
        int Ny = 0;
        int pixel_size = 0;

        attr->read(*type, &Nx);

        ccd->openAttribute("Ny").read(*type, &Ny);
        ccd->openAttribute("pixelsize").read(*type, &pixel_size);

        this->ccd = new CCD(Nx, Ny, pixel_size);
        delete ccd;
    }


    delete type;
    delete attr;

    delete h5file;

    if (file_contains_psfs) {
        this->psfs = new PSF_ZEMAX(path, fiber_number);
    }
}

double MatrixSimulator::get_maximum_wavelength() {
    return this->wavelength_limit_max;
}

double MatrixSimulator::get_minimum_wavelength() {
    return this->wavelength_limit_min;
}

double MatrixSimulator::get_blaze() {
    return this->spec_info.blaze;
}

double MatrixSimulator::get_gpmm() {
    return this->spec_info.gpmm;
}

void MatrixSimulator::set_telescope(Telescope *telescope) {
    this->telescope = *telescope;
}

void MatrixSimulator::calc_splines() {
    this->tr_p.clear();
    this->tr_q.clear();
    this->tr_phi.clear();
    this->tr_r.clear();
    this->tr_tx.clear();
    this->tr_ty.clear();

    for (int o = 0; o < this->orders.size(); ++o) {
        std::vector<double> wl;
        std::vector<double> p;
        std::vector<double> q;
        std::vector<double> r;
        std::vector<double> phi;
        std::vector<double> tx;
        std::vector<double> ty;

        for (auto rt: this->raw_transformations[o + this->min_order]) {
            wl.push_back(rt.wavelength);
            p.push_back(rt.decomposed_matrix[0]);
            q.push_back(rt.decomposed_matrix[1]);
            r.push_back(rt.decomposed_matrix[2]);
            phi.push_back(rt.decomposed_matrix[3]);
            tx.push_back(rt.decomposed_matrix[4]);
            ty.push_back(rt.decomposed_matrix[5]);
        }
        tk::spline s_p;
        s_p.set_points(wl, p);
        this->tr_p.push_back(s_p);

        tk::spline s_q;
        s_q.set_points(wl, q);
        this->tr_q.push_back(s_q);

        tk::spline s_r;
        s_r.set_points(wl, r);
        this->tr_r.push_back(s_r);

        tk::spline s_phi;
        s_phi.set_points(wl, phi);
        this->tr_phi.push_back(s_phi);

        tk::spline s_tx;
        s_tx.set_points(wl, tx);
        this->tr_tx.push_back(s_tx);

        tk::spline s_ty;
        s_ty.set_points(wl, ty);
        this->tr_ty.push_back(s_ty);

    }

}

void MatrixSimulator::set_wavelength(int N) {
    this->sim_wavelength.clear();
    for (int o = 0; o < this->orders.size(); ++o)
        this->sim_wavelength.push_back(std::vector<double>(N));

#pragma omp parallel for
    for (int o = 0; o < this->orders.size(); ++o) {
        double min_wl = this->raw_transformations[o + this->min_order].front().wavelength;
        double max_wl = this->raw_transformations[o + this->min_order].back().wavelength;
        double dl = (max_wl - min_wl) / N;
        for (int i = 0; i < N; ++i) {
            this->sim_wavelength[o][i] = min_wl + dl * i;
        }
    }
}

void MatrixSimulator::set_wavelength(std::vector<double> wavelength) {
    this->sim_wavelength.clear();
    for (int o = 0; o < this->orders.size(); ++o) {
        double min_wl = this->raw_transformations[o + this->min_order].front().wavelength;
        double max_wl = this->raw_transformations[o + this->min_order].back().wavelength;
        std::vector<double> wl_in_order;

        for (auto wl: wavelength) {
            if ((wl > min_wl) && (wl < max_wl)) {
                wl_in_order.push_back(wl);
            }
        }
        this->sim_wavelength.push_back(wl_in_order);
    }
}


std::array<float, 6> MatrixSimulator::get_transformation_matrix_lookup(int o, double wavelength) {
    std::vector<float> parameters;
    int idx_matrix = floor((wavelength - this->sim_wavelength[o].front()) / this->sim_matrix_dwavelength[o]);

    parameters.push_back(this->sim_p[o][idx_matrix]);
    parameters.push_back(this->sim_q[o][idx_matrix]);
    parameters.push_back(this->sim_r[o][idx_matrix]);
    parameters.push_back(this->sim_phi[o][idx_matrix]);
    parameters.push_back(this->tr_tx[o](wavelength));
    parameters.push_back(this->tr_ty[o](wavelength));
    return compose_matrix(parameters);

}

std::array<float, 6> MatrixSimulator::get_transformation_matrix(int o, double wavelength) {
    std::vector<float> parameters;
    parameters.push_back(this->tr_p[o](wavelength));
    parameters.push_back(this->tr_q[o](wavelength));
    parameters.push_back(this->tr_r[o](wavelength));
    parameters.push_back(this->tr_phi[o](wavelength));
    parameters.push_back(this->tr_tx[o](wavelength));
    parameters.push_back(this->tr_ty[o](wavelength));
    return compose_matrix(parameters);
}

int MatrixSimulator::simulate(double t, unsigned long seed) {
    this->set_efficiencies(this->efficiencies);

    this->prepare_sources(this->sources); //put area eventually in sources

    this->prepare_psfs(1000);

    this->prepare_matrix_lookup(1000);

    std::vector<Histogram> wl_s(orders.size());
    std::vector<int> N_photons(orders.size());

    double psf_scaling = (*this->ccd->get_pixelsize() / this->psfs->pixelsampling);

    std::cout << "Number of photons per order:" << std::endl;
    int N_tot = 0;
    for (int o = 0; o < this->orders.size(); ++o) {
        if (!sim_wavelength[o].empty()) {
            std::vector<double> a(sim_wavelength[o].begin(), sim_wavelength[o].end()); //units are um
            std::vector<double> b(sim_spectra_time_efficieny[o].begin(),
                                  sim_spectra_time_efficieny[o].end()); //units are uW per um

            wl_s[o] = Histogram(a, b);
            wl_s[o].mode = this->mode;
            //units are assumed to be t=[s], area=[m^2], wl_s.dflux=[Num of Photons]/([s] * [m^2] * [um]), wl_s.Calc_flux = [Num of Photons]/([s]*[m^2])
            N_photons[o] = static_cast<int>(floor(wl_s[o].Calc_flux() * t));
            N_tot += N_photons[o];
            //this->telescope->get_area();
            //t*area*wl_s[o].Calc_flux()
            std::cout << "Order " << o + this->min_order << ": " << N_photons[o] << std::endl;
        } else {
            N_photons[o] = 0;
            std::cout << "Order " << o + this->min_order << ": " << N_photons[o] << std::endl;
        }
    }
    std::cout << "Total number of photons:" << N_tot << std::endl;
    std::cout << "Start tracing ..." << std::endl;
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    std::vector<int> local_data = this->ccd->data;
    for (int o = 0; o < this->orders.size(); ++o) {
        cout << "Simulating Order... " << o + this->min_order << "/" << this->max_order << std::endl;
#pragma omp parallel
        {
            std::uniform_real_distribution<float> dist(0.0, 1.0);
            std::discrete_distribution<int> disf(wl_s[o].intensity.begin(), wl_s[o].intensity.end());
//            std::uniform_real_distribution<float> rgx(0., this->slit->slit_sampling);
//            std::uniform_real_distribution<float> rgy(0., this->slit->slit_sampling * this->slit->h / this->slit->w);
            std::uniform_real_distribution<float> rgx(0., 1.);
            std::uniform_real_distribution<float> rgy(0., 1.);
//            std::default_random_engine gen;
            std::mt19937 gen;
            if (seed == 0) {
                std::random_device rd;
                gen.seed(rd());
            } else {
#if defined(_OPENMP)
                //TODO: better handling of seeds for parallel execution.. (i.e. use random seed sequence)
                gen.seed(seed + omp_get_thread_num());
#else
                gen.seed(seed);
#endif
            }

#pragma omp for reduction(vec_int_plus : local_data)
            for (int i = 0; i < N_photons[o]; ++i) {
                //double wl = wl_s[o].Sample(dist(gen));
                double wl = wl_s[o].event[disf(gen)];
                //cout<<wl<<endl;

                int idx_matrix = (floor((wl - this->sim_wavelength[o].front()) / this->sim_matrix_dwavelength[o]));

                float x = rgx(gen);
                float y = rgy(gen);

                float newx = (sim_m00[o][idx_matrix] * x + sim_m01[o][idx_matrix] * y + this->tr_tx[o](wl));
                float newy = (sim_m10[o][idx_matrix] * x + sim_m11[o][idx_matrix] * y + this->tr_ty[o](wl));

                //            this is not accurate enough.  probably fixable by linear interpolation of tx and ty
                //            float newx = (sim_m00[o][idx_matrix] * x + sim_m01[o][idx_matrix] * y + sim_tx[o][idx_matrix]);
                //            float newy = (sim_m10[o][idx_matrix]  * x + sim_m11[o][idx_matrix]  * y + sim_ty[o][idx_matrix]);

                // lookup psf
                auto idx_psf = static_cast<int> (floor(
                        (wl - this->sim_wavelength[o].front()) / this->sim_psfs_dwavelength[o]));

                std::uniform_real_distribution<float> abr_x(0, this->sim_psfs[o][idx_psf].cols - 1);
                std::uniform_real_distribution<float> abr_y(0, this->sim_psfs[o][idx_psf].rows - 1);
                std::uniform_real_distribution<float> abr_z(0, 1);
                float abx = abr_x(gen);
                float aby = abr_y(gen);
                float abz = abr_z(gen);
                while (abz > this->sim_psfs[o][idx_psf].data[floor(aby)][floor(abx)]) {
                    abx = abr_x(gen);
                    aby = abr_y(gen);
                    abz = abr_z(gen);
                }

                newx += (abx - this->sim_psfs[o][idx_psf].cols / 2.) / psf_scaling;
                newy += (aby - this->sim_psfs[o][idx_psf].rows / 2.) / psf_scaling;

                if (newx > 0 && newx < this->ccd->Nx && newy > 0 && newy < this->ccd->Ny)
                    local_data[(int) floor(newx) + this->ccd->Nx * (int) floor(newy)] += 1;
            }
        };
        this->ccd->data = local_data;
    };
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000000.;
    std::cout << "Duration Tracing: \t" << duration << " s" << std::endl;
}

void MatrixSimulator::set_efficiencies(std::vector<Efficiency *> &efficiencies) {
    std::cout << std::endl << "Set Efficiencies ... " << std::endl;
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    for (int o = 0; o < this->orders.size(); ++o) {
        this->sim_efficiencies.push_back(std::vector<double>(this->sim_wavelength[o].size(), 1.));
        this->sim_total_efficiency_per_order.push_back(0.);
    }

#pragma omp parallel for
    for (int o = 0; o < this->orders.size(); ++o) {
        for (auto &e : efficiencies) {
            std::vector<double> sim_eff = e->get_efficiency(o + this->min_order, this->sim_wavelength[o]);
            for (int i = 0; i < this->sim_efficiencies[o].size(); ++i) {
                this->sim_efficiencies[o][i] *= sim_eff[i];
                this->sim_total_efficiency_per_order[o] += sim_eff[i];
            }
        }
    }
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000000.;
    std::cout << "Duration: \t" << duration << " s" << std::endl;
}

void MatrixSimulator::prepare_psfs(int N) {

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Prepare PSFs ..." << std::endl;

    for (int o = 0; o < this->orders.size(); ++o) {
        this->sim_psfs.emplace_back(std::vector<Matrix>(N));
        this->sim_psfs_wavelength.emplace_back(std::vector<double>(N));
        this->sim_psfs_dwavelength.push_back(0.);
    }

#pragma omp parallel for
    for (int o = 0; o < this->orders.size(); ++o) {
        double min_wl = this->raw_transformations[o + this->min_order].front().wavelength;
        double max_wl = this->raw_transformations[o + this->min_order].back().wavelength;
        double dl = (max_wl - min_wl) / N;
        this->sim_psfs_dwavelength[o] = dl;
        for (int i = 0; i < this->sim_psfs_wavelength[o].size(); ++i) {
            this->sim_psfs_wavelength[o][i] = min_wl + i * dl;
            this->sim_psfs[o][i] = Matrix(this->psfs->get_PSF(o + this->min_order, min_wl + i * dl));
        }
    }
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000000.;
    std::cout << "Duration: \t" << duration << " s" << std::endl;
}


void MatrixSimulator::prepare_matrix_lookup(int N) {

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Prepare matrix elements ..." << std::endl;

    for (int o = 0; o < this->orders.size(); ++o) {
        this->sim_matrix_wavelength.push_back(std::vector<double>(N));
        this->sim_matrix_dwavelength.push_back(0.);
    }

    for (int o = 0; o < this->orders.size(); ++o) {
        double min_wl = this->raw_transformations[o + this->min_order].front().wavelength;
        double max_wl = this->raw_transformations[o + this->min_order].back().wavelength;
        double dl = (max_wl - min_wl) / N;
        this->sim_matrix_dwavelength[o] = dl;
        std::vector<double> p;
        std::vector<double> q;
        std::vector<double> r;
        std::vector<double> phi;
        std::vector<double> tx;
        std::vector<double> ty;

        std::vector<double> m00;
        std::vector<double> m10;
        std::vector<double> m01;
        std::vector<double> m11;


        for (int i = 0; i < this->sim_matrix_wavelength[o].size(); ++i) {
            this->sim_matrix_wavelength[o][i] = min_wl + i * dl;

            p.push_back(this->tr_p[o](this->sim_matrix_wavelength[o][i]));
            q.push_back(this->tr_q[o](this->sim_matrix_wavelength[o][i]));
            r.push_back(this->tr_r[o](this->sim_matrix_wavelength[o][i]));
            phi.push_back(this->tr_phi[o](this->sim_matrix_wavelength[o][i]));
            tx.push_back(this->tr_tx[o](this->sim_matrix_wavelength[o][i]));
            ty.push_back(this->tr_ty[o](this->sim_matrix_wavelength[o][i]));

            std::vector<float> params;
            params.push_back(p[i]);
            params.push_back(q[i]);
            params.push_back(r[i]);
            params.push_back(phi[i]);
            params.push_back(tx[i]);
            params.push_back(ty[i]);
            std::array<float, 6> tm = compose_matrix(params);

            m00.push_back(tm[0]);
            m01.push_back(tm[1]);
            m10.push_back(tm[3]);
            m11.push_back(tm[4]);
        }
        this->sim_p.push_back(p);
        this->sim_q.push_back(q);
        this->sim_r.push_back(r);
        this->sim_phi.push_back(phi);
        this->sim_tx.push_back(tx);
        this->sim_ty.push_back(ty);
        this->sim_m00.push_back(m00);
        this->sim_m01.push_back(m01);
        this->sim_m10.push_back(m10);
        this->sim_m11.push_back(m11);


    }
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000000.;
    std::cout << "Duration: \t" << duration << " s" << std::endl;
}

void MatrixSimulator::prepare_sources(std::vector<Source *> sources) {
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Prepare sources ..." << std::endl;
    this->sim_spectra.clear();
    // allocate memory so we can use omp parallel for filling
    for (int o = 0; o < this->orders.size(); ++o) {
        this->sim_spectra.push_back(std::vector<float>(this->sim_wavelength[o].size()));
        this->sim_spectra_time_efficieny.push_back(std::vector<float>(this->sim_wavelength[o].size()));
    }

#pragma omp parallel for
    for (int o = 0; o < this->orders.size(); ++o) {
        {
            for (auto &s : sources) {
                std::vector<float> spectrum = s->get_spectrum(this->sim_wavelength[o]);
                for (int i = 0; i < spectrum.size(); ++i) {
                    this->sim_spectra[o][i] = spectrum[i] * telescope.get_area();
                    this->sim_spectra_time_efficieny[o][i] = spectrum[i] * this->sim_efficiencies[o][i];
                }
            }
        }
    }
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000000.;
    std::cout << "Duration: \t" << duration << " s" << std::endl;
}

void MatrixSimulator::add_efficiency(Efficiency *eff) {
    this->efficiencies.push_back(eff);
}


void MatrixSimulator::set_ccd(CCD *ccd) {
    this->ccd = ccd;
}

void MatrixSimulator::set_slit(Slit *slit) {
    this->slit = slit;
}

void MatrixSimulator::set_psfs(PSF *psfs) {
    this->psfs = psfs;
}

void MatrixSimulator::set_source(Source *src) {
    this->sources.push_back(src);
}

void MatrixSimulator::save_to_hdf(const std::string filename, bool downsample, bool bleed, bool overwrite) {
    this->ccd->save_to_hdf(filename, downsample, bleed, overwrite);
}

void MatrixSimulator::save_to_fits(const std::string filename, bool downsample, bool bleed, bool overwrite) {
    this->ccd->save_to_fits(filename, overwrite);
}

void MatrixSimulator::save_1d_to_fits(const std::string filename) {
    CCfits::FITS infile(filename.c_str(), CCfits::Write);
    string hduName("1Dspectrum_Fiber_" + std::to_string(this->fiber_number));
    int n_orders = this->orders.size();

    CCfits::Table *newTable = infile.addTable(hduName, 2);

    for (auto const &o: this->orders) {
        newTable->addColumn(CCfits::Tdouble, "Order_" + std::to_string(o), this->sim_spectra[o].size(), "");
        newTable->column("Order_" + std::to_string(o)).write(this->sim_spectra[o], 1, 2);
        newTable->column("Order_" + std::to_string(o)).write(this->sim_wavelength[o], 1, 1);
    }

}

void MatrixSimulator::add_background(double bias, double noise, unsigned long seed) {
    std::mt19937 gen;

    if (seed == 0) {
        std::random_device rd;
        gen.seed(rd());
    } else {
        gen.seed(seed);
    }
    std::normal_distribution<double> dis(bias, noise);

    for (int i = 0; i < this->ccd->Nx * this->ccd->Ny; ++i) {
        this->ccd->data[i] += dis(gen);
    }
}