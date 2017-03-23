#include "matrixsimulator.h"
#include <iostream>
#include <fstream>
#include <iterator>
#include "spline.h"
// #include <string>
// #include <stdlib.h>
#include "H5Cpp.h"
#include <math.h>
// using namespace std
#include "helper.h"
#include <opencv2/imgproc.hpp>
#include <hdf5_hl.h>
#include <CCfits/CCfits>
#include <CCfits/FITS.h>
//#include <curlpp/cURLpp.hpp>
//#include <curlpp/Easy.hpp>
//#include <curlpp/Options.hpp>
#include <string>
#include <fstream>

#ifdef USE_GPU
#include "opencv2/gpu/gpu.hpp"
#endif

void print_transformation_matrix(cv::Mat tm){
    std::cout << tm.at<double>(0,0) << "\t" << tm.at<double>(0,1) << "\t" << tm.at<double>(0,2) << "\n";
    std::cout << tm.at<double>(1,0) << "\t" << tm.at<double>(1,1) << "\t" << tm.at<double>(1,2) << std::endl;

}

typedef struct transformation_hdf
{
    float rotation;
    float scale_x;
    float scale_y;
    float shear;
    float translation_x;
    float translation_y;
    float wavelength;
} transformation_hdf;

MatrixSimulator::MatrixSimulator()
{
}

//int download_spectrograph_model(){
//    try
//    {
//        std::string url = "http://www.example.com/";
//        std::ofstream out("output");
//
//        curlpp::Cleanup cleanup;
//        curlpp::Easy request;
//
//        request.setOpt(new curlpp::options::Url(url));
//        request.setOpt(new curlpp::options::WriteStream(&out));
//
//        request.perform();
//    }
//    catch( cURLpp::RuntimeError &e )
//    {
//        std::cerr << e.what() << std::endl;
//        return 1;
//    }
//    catch( cURLpp::LogicError &e )
//    {
//        std::cerr << e.what() << std::endl;
//        return 1;
//    }
//}

void MatrixSimulator::load_spectrograph_model(std::string path, int fiber_number, bool keep_ccd) {
    if (path.find('.hdf') != std::string::npos) {
        std::cout << "found!" << '\n';
    } else
    {
//        download_spectrograph_model()https://github.com/Stuermer/EchelleSimulator/blob/master/data/MaroonX.hdf
    }
    const H5std_string filename(path);
    this->efficiencies.clear();
    this->orders.clear();
    this->raw_transformations.clear();
    this->sim_wavelength.clear();
    this->sim_spectra.clear();
    this->sim_efficiencies.clear();
    this->sim_matrices.clear();

    // open file readonly
    H5::H5File *h5file = new H5::H5File(filename, H5F_ACC_RDONLY);

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
    int oversampling = 0;
    int sampling_input_x;
    int number_of_points;

    attr->read(*type, &field_height);
    fiber->openAttribute("field_with").read(*type, &field_width);
    attr = new H5::Attribute(fiber->openAttribute("oversampling_output"));
    type = new H5::DataType(attr->getDataType());
    attr->read(*type, &oversampling);
    fiber->openAttribute("sampling_input_x").read(*type, &sampling_input_x);
    fiber->openAttribute("MatricesPerOrder").read(*type, &number_of_points);

    this->slit = new Slit(field_width, field_height, sampling_input_x);

    std::vector<std::string> group_names;
    herr_t idx = H5Literate(fiber->getId(), H5_INDEX_NAME, H5_ITER_INC, NULL, file_info, &group_names);
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

    for (auto gn: group_names) {
        if (gn.find("psf") == std::string::npos) {
            transformation_hdf data[number_of_points];
            std::string tablename = "fiber_" + std::to_string(fiber_number) + "/" + gn;
            int order = std::stoi(gn.substr(5));
            herr_t read = H5TBread_table(h5file->getId(), tablename.c_str(), dst_size, dst_offset, dst_sizes, &data);
//            std::cout << order << "\t" << read << "\t" << data[0].translation_x << std::endl;
            for (int i = 0; i < number_of_points; ++i) {
                raw_transformation t;
                t.wavelength = data[i].wavelength;
                std::vector<double> result;
                // return <sx, sy, shear, rot, tx ,ty>
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

    for (auto imap: this->raw_transformations)
        this->orders.push_back(imap.first);
    this->calc_splines();

    if (!keep_ccd) {
        // read in CCD information
        H5::Group *ccd = new H5::Group(h5file->openGroup("CCD"));
        attr = new H5::Attribute(ccd->openAttribute("Nx"));
        type = new H5::DataType(attr->getDataType());
        int Nx = 0;
        int Ny = 0;
        int pixelsize = 0;

        attr->read(*type, &Nx);

        ccd->openAttribute("Ny").read(*type, &Ny);
        ccd->openAttribute("pixelsize").read(*type, &pixelsize);

        this->ccd = new CCD(Nx, Ny, 1, this->slit->slit_image.type());
        delete ccd;
    }


    delete type;
    delete attr;

    delete h5file;

    if (file_contains_psfs){
            this->psfs = new PSF_ZEMAX(path, fiber_number);
        }
}

void MatrixSimulator::read_transformations(std::string path)
{
    std::ifstream       file(path.c_str());

    for(CSVIterator loop(file); loop != CSVIterator(); ++loop)
    {
        raw_transformation t;
        int order = std::stoi((*loop)[0]);

        t.order = std::stoi((*loop)[0]);
        t.wavelength = std::stod(std::string((*loop)[1]));

        t.transformation_matrix(0,0) = std::stod((*loop)[2]);
        t.transformation_matrix(0,1) = std::stod((*loop)[3]);
        t.transformation_matrix(0,2)= std::stod((*loop)[4]);

        t.transformation_matrix(1,0) = std::stod((*loop)[5]);
        t.transformation_matrix(1,1) = std::stod((*loop)[6]);
        t.transformation_matrix(1,2) = std::stod((*loop)[7]);

        // std::cout << t.order << " " << t.wavelength << " " << t.transformation_matrix.at<double>(0,0) <<"\n";
        t.decomposed_matrix = decompose_matrix(t.transformation_matrix);
        // std::vector<double> tmp = t.decomposed_matrix;
        // std::cout << tmp[0]<< " " << tmp[1]<< " " << tmp[2]<< " " << tmp[3];
        this->raw_transformations[order].push_back(t);
        // this->raw_transformations.push_back(t);

        // std::cout << std::to_string(t.order) << "  " << t.wavelength << "  "<< t.transformation_matrix.at<double>(0,0) << "\n";
    }

    for(auto imap: this->raw_transformations)
        this->orders.push_back(imap.first);
    this->calc_splines();
}


double MatrixSimulator::get_blaze() {
    return this->spec_info.blaze;
}

double MatrixSimulator::get_gpmm() {
    return this->spec_info.gpmm;
}

void MatrixSimulator::calc_splines(){
    for(auto o: this->orders)
    {
        std::vector<double> wl;
        std::vector<double> p;
        std::vector<double> q;
        std::vector<double> r;
        std::vector<double> phi;
        std::vector<double> tx;
        std::vector<double> ty;

        for(auto rt: this->raw_transformations[o])
        {
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
        this->tr_p[o] = s_p;

        tk::spline s_q;
        s_q.set_points(wl, q);
        this->tr_q[o] = s_q;

        tk::spline s_r;
        s_r.set_points(wl, r);
        this->tr_r[o] = s_r;

        tk::spline s_phi;
        s_phi.set_points(wl, phi);
        this->tr_phi[o] = s_phi;

        tk::spline s_tx;
        s_tx.set_points(wl, tx);
        this->tr_tx[o] = s_tx;

        tk::spline s_ty;
        s_ty.set_points(wl, ty);
        this->tr_ty[o] = s_ty;

    }

}

void MatrixSimulator::set_wavelength(int N){
    this->sim_wavelength.clear();

//#pragma omp parallel for
    for(auto o: this->orders)
    {
        double min_wl = this->raw_transformations[o].front().wavelength;
        double max_wl = this->raw_transformations[o].back().wavelength;
        double dl = (max_wl - min_wl) / N;
        for(int i=0; i<N; ++i)
        {
            this->sim_wavelength[o].push_back(min_wl+dl*i);
        }
    }
    this->calc_sim_matrices();
}

void MatrixSimulator::set_wavelength(std::vector<double> wavelength){
    this->sim_wavelength.clear();
    for(auto o: this->orders)
    {
        double min_wl = this->raw_transformations[o].front().wavelength;
        double max_wl = this->raw_transformations[o].back().wavelength;
        std::vector<double> wl_in_order;

        for (auto wl: wavelength)
        {
            if ((wl>min_wl) && (wl<max_wl))
            {
                wl_in_order.push_back(wl);
            }
        }
        this->sim_wavelength.insert(std::pair<int, std::vector<double>> (o, wl_in_order));
    }
    this->calc_sim_matrices();
}

void MatrixSimulator::calc_sim_matrices(){
    std::cout << "Calculate Matrices ..." << std::endl;
    this->sim_matrices.clear();
#pragma omp parallel for
    for(int o=this->orders.front(); o<this->orders.back()+1; ++o)
    {
//        std::cout<< o << std::endl;
        for(auto w: this->sim_wavelength[o])
        {
            this->sim_matrices[o].push_back(this->get_transformation_matrix(o, w));
        }
        /*
        std::cout <<  this->sim_matrices[o][0].at<double>(0,0)<<"\t";
        std::cout <<  this->sim_matrices[o][0].at<double>(0,1)<<"\t";
        std::cout <<  this->sim_matrices[o][0].at<double>(0,2)<<"\n";
        std::cout <<  this->sim_matrices[o][0].at<double>(1,0)<<"\t";
        std::cout <<  this->sim_matrices[o][0].at<double>(1,1)<<"\t";
        std::cout <<  this->sim_matrices[o][0].at<double>(1,2)<<"\n\n";
        */
    }
}

Matrix23f MatrixSimulator::get_transformation_matrix(int o, double wavelength){
    std::vector<double> parameters;
    parameters.push_back(this->tr_p[o](wavelength));
    parameters.push_back(this->tr_q[o](wavelength));
    parameters.push_back(this->tr_r[o](wavelength));
    parameters.push_back(this->tr_phi[o](wavelength));
    parameters.push_back(this->tr_tx[o](wavelength));
    parameters.push_back(this->tr_ty[o](wavelength));
    return compose_matrix(parameters);
}

#ifdef USE_GPU
cv::gpu::GpuMat MatrixSimulator::transform_slit(cv::gpu::GpuMat& slit_image, cv::Mat& transformation_matrix, double weight){
    int n_rows = abs(round(slit_image.rows * 1.5 * transformation_matrix.at<double>(1,1)));
    int n_cols = abs(round(slit_image.cols * 4. * transformation_matrix.at<double>(0,0)));

    cv::Mat tm = transformation_matrix.clone();
    int tx_int = floor(tm.at<double>(0,2));
    int ty_int = floor(tm.at<double>(1,2));
    tm.at<double>(0,2) += -tx_int + n_cols/2.;
    tm.at<double>(1,2) += -ty_int + n_rows * 0.75;

    // cv::Mat warp_dst_cpu = cv::Mat::zeros( n_rows, n_cols, slit_image.type() );
    cv::gpu::GpuMat warp_dst=cv::gpu::GpuMat(n_rows, n_cols, slit_image.type() );
    //warp_dst.upload(warp_dst_cpu);

    cv::gpu::warpAffine(slit_image.clone(), warp_dst, tm, warp_dst.size() );
    return warp_dst;
}
#endif
double tmp_image_height(int rows, double sy)
{
    return abs(rows * sy)+16.;
}

double tmp_image_width(int cols, double sx)
{
    return abs(cols * sx)+16.;
}

cv::Mat MatrixSimulator::transform_slit(cv::Mat& slit_image, cv::Mat& transformation_matrix, double weight){
    double a,b,c,d;
    a= transformation_matrix.at<double>(0,0);
    b= transformation_matrix.at<double>(1,0);
    c= transformation_matrix.at<double>(0,1);
    d= transformation_matrix.at<double>(1,1);

    double x00 = transformation_matrix.at<double>(0,2);
    double y00 = transformation_matrix.at<double>(1,2);
    double x11 = a*this->slit->w_px+b*this->slit->h_px+x00;
    double y11 = c*this->slit->w_px+d*this->slit->h_px+y00;

    int xlu = std::min(floor(x00),floor(x11));
    int ylu = std::min(floor(y00),floor(y11));
    int xro = std::max(floor(x00),floor(x11));
    int yro = std::max(floor(y00),floor(y11));

    int n_rows = yro - ylu + 16;
    int n_cols = xro - xlu + 16;

    cv::Mat tm = transformation_matrix.clone();

    int offset_x = xlu - 8;
    int offset_y = ylu - 8;


    tm.at<double>(0,2) -= offset_x;
    tm.at<double>(1,2) -= offset_y;

    double det = cv::determinant(transformation_matrix.colRange(0,2));
    cv::Mat warp_dst;
    weight /= det;
    warp_dst = cv::Mat::zeros( n_rows, n_cols, slit_image.type() );
    cv::warpAffine(slit_image.clone()*weight, warp_dst, tm, warp_dst.size() );

    return warp_dst;
}

#ifdef USE_GPU
int MatrixSimulator::simulate_order(int order, cv::gpu::GpuMat& slit_image, cv::gpu::GpuMat& output_image, bool aberrations)
{
    // cv::Mat img = cv::Mat::zeros(4096*3, 4096*3, CV_64FC1);
    int n = this->sim_matrices[order].size();


    for(int i=0; i<this->sim_matrices[order].size(); ++i){
        cv::Mat tr=this->sim_matrices[order][i];
        double weight = this->sim_efficiencies[order][i];
        double wavelength = this->sim_wavelength[order][i];
        cv::gpu::GpuMat tmp = cv::gpu::GpuMat();
        cv::Mat tmp2;
        tmp = transform_slit(slit_image, tr, this->sim_efficiencies[order][i]);

        if (aberrations){
            PSF * psf_p = dynamic_cast<PSF*> (this->psfs);
            cv::Mat psf = psf_p->get_PSF(order, wavelength);
            cv::flip(psf, psf, -1);
//          cv::gpu::GpuMat psf_gpu =  cv::gpu::GpuMat(psf);
            cv::gpu::filter2D(tmp, tmp, -1, psf, cvPoint(psf.cols/2, psf.rows/2));
//          cv::filter2D(tmp, tmp, -1, psf, cvPoint(psf.cols/2, psf.rows/2));
        }

        int tx_int = floor(tr.at<double>(0,2));
        int ty_int = floor(tr.at<double>(1,2));

        int n_rows = abs(round(slit_image.rows * 1.5 * tr.at<double>(1,1)));
        int n_cols = abs(round(slit_image.cols * 4. * tr.at<double>(0,0)));


        int left_margin = std::max(0, tx_int);
        int right_margin = std::min(output_image.cols, tx_int+n_cols);
        int width = right_margin - left_margin;

        cv::gpu::multiply(tmp, weight, tmp);

        if (width>0 && width>=n_cols && (ty_int>0) && (ty_int+n_rows<output_image.rows) && tx_int>0 && tx_int+n_cols<output_image.cols)
        {
            cv::gpu::GpuMat roi = output_image.rowRange(ty_int, ty_int + n_rows).colRange(tx_int, tx_int+n_cols);
            cv::gpu::add(roi, tmp, roi );
        }
    }
    return 0;
}
#endif

int MatrixSimulator::photon_order(int N_photons) {
    std::cout <<"Start tracing ..." <<std::endl;


    // global img
    std::vector<uint16_t> img(this->ccd->data.rows * this->ccd->data.cols, 0);

//    #pragma omp parallel
//    {
//         private img
//        std::vector<uint16_t> img_private(this->ccd->data.rows * this->ccd->data.cols, 0);
        #pragma omp parallel for
        for(int o=this->orders.front(); o<this->orders.back(); ++o) {
            std::random_device rd;
            std::default_random_engine gen(rd());
            std::uniform_real_distribution<> dis_slitx(0., this->slit->slit_sampling);
            std::uniform_real_distribution<> dis_slity(0., this->slit->slit_sampling * this->slit->h/this->slit->w);
            std::cout << "Order... " << o << std::endl;
            double min_wl = this->raw_transformations[o].front().wavelength;
            double max_wl = this->raw_transformations[o].back().wavelength;

            int N = 10000;
            double dl = (max_wl - min_wl) / N;
            std::vector<double> wl;
            for(int i=0; i<N; ++i)
            {
                wl.push_back(min_wl+dl*i);
            }

            std::vector<double> specD;
            std::vector<float> spec;
            specD = this->sources[0]->get_spectrum(wl);

            for(int j=0; j<specD.size(); ++j)
                spec.push_back(float(specD[j]));

            std::vector<double> totaleffD(wl.size(), 1.);
            std::vector<float> eff;
            for(auto& e : this->efficiencies)
            {
                std::vector<double> sim_eff = e->get_efficieny(o, wl);
                for(int i=0; i<totaleffD.size(); ++i)
                {
                    totaleffD[i] *= sim_eff[i];
                }
            }

            for(int k=0; k<totaleffD.size(); ++k)
                eff.push_back(float(totaleffD[k]) * float(specD[k]));

            std::piecewise_linear_distribution<> dis(wl.begin(), wl.end(), eff.begin());

//            std::uniform_real_distribution<> dis(min_wl, max_wl);

            for (int i = 0; i < N_photons; ++i) {
                double wl = dis(gen);
                Matrix23f tm = this->get_transformation_matrix(o, wl);


                float x = dis_slitx(gen);
                float y = dis_slity(gen);

                float newx = (tm(0, 0) * x + tm(0, 1) * y + tm(0, 2)) / 3.;
                float newy = (tm(1, 0) * x + tm(1, 1) * y + tm(1, 2)) / 3.;

                cv::Mat psf = this->psfs->get_PSF(o, wl);
                std::uniform_int_distribution<> abr_x(0, psf.cols);
                std::uniform_int_distribution<> abr_y(0, psf.rows);
                double mx;
                cv::minMaxIdx(psf,NULL, &mx);
                std::uniform_real_distribution<> abr_z(0,mx);
                float abx=0.;
                float aby=0.;
                float abz=1.;
                while(abz > psf.at<double>(floor(aby), floor(abx))) {
                    abx = abr_x(gen);
                    aby = abr_y(gen);
                    abz = abr_z(gen);
                }
//                std::cout<<(abx-psf.cols/2.)/3.<<"\t"<<(aby-psf.rows/2.)/3.<< std::endl;
                newx += (abx-psf.cols/2.)/3.;
                newy += (aby-psf.rows/2.)/3.;

                if (newx > 0 && newx < this->ccd->data.cols && newy > 0 && newy < this->ccd->data.rows)
                    img[floor(newx) + floor(newy) * this->ccd->data.cols] += 1;
                //                this->ccd->data.at<double>(floor(newy), floor(newx)) += 1.;

            }
        }
//        #pragma omp critical
//        {
//            for(int i=0; i<img_private.size(); ++i)
//                img[i] += img_private[i];
//        }

//    }

    for(int xx=0; xx<this->ccd->data.cols; ++xx)
        for(int yy=0; yy<this->ccd->data.rows; ++yy)
            this->ccd->data.at<double>(xx,yy) += img[yy+xx*this->ccd->data.cols];

}

int MatrixSimulator::simulate_order(int order, cv::Mat& slit_image, cv::Mat& output_image, bool aberrations)
{

    // cv::Mat img = cv::Mat::zeros(4096*3, 4096*3, CV_64FC1);
    int n = this->sim_matrices[order].size();
    std::cout << "Simulate order " << order << " number of wavelength " << n << std::endl;

    for(int i=0; i<n; ++i)
    {
        Matrix23f tr_eigen = this->sim_matrices[order][i];
        cv::Mat tr = cv::Mat(2,3, CV_64FC1);
        tr.at<double>(0,0) = tr_eigen(0,0);
        tr.at<double>(0,1) = tr_eigen(0,1);
        tr.at<double>(0,2) = tr_eigen(0,2);
        tr.at<double>(1,0) = tr_eigen(1,0);
        tr.at<double>(1,1) = tr_eigen(1,1);
        tr.at<double>(1,2) = tr_eigen(1,2);


        double weight = this->sim_efficiencies[order][i];
        double wavelength = this->sim_wavelength[order][i];
//        weight *= this->sim_spectra[order][i];
        cv::Mat tmp = transform_slit(slit_image, tr, weight);

        if (aberrations){
            PSF * psf_p = dynamic_cast<PSF *> (this->psfs);
            cv::Mat psf = psf_p->get_PSF(order, wavelength);
            cv::flip(psf, psf, -1);
            cv::filter2D(tmp, tmp, -1, psf, cvPoint(psf.cols/2, psf.rows/2));
        }


        double a,b,c,d;
        a = tr_eigen(0,0);
        b = tr_eigen(1,0);
        c = tr_eigen(0,1);
        d = tr_eigen(1,1);

        double x00 = tr_eigen(0,2);
        double y00 = tr_eigen(1,2);
        double x11 = a*this->slit->w_px+b*this->slit->h_px+x00;
        double y11 = c*this->slit->w_px+d*this->slit->h_px+y00;

        int xlu = std::min(floor(x00),floor(x11));
        int ylu = std::min(floor(y00),floor(y11));
        int xro = std::max(floor(x00),floor(x11));
        int yro = std::max(floor(y00),floor(y11));

        int n_rows = yro - ylu + 16;
        int n_cols = xro - xlu + 16;

//        cv::Mat tm = tr.clone();

        int tx_int = xlu - 8;
        int ty_int = ylu - 8;


        int left_margin = std::max(0, tx_int);
        int right_margin = std::min(output_image.cols, tx_int+n_cols);
        int width = right_margin - left_margin;

        if (width>0 && width>=n_cols && (ty_int>0) && (ty_int+n_rows<output_image.rows) && tx_int>0 && tx_int+n_cols<output_image.cols)
        {
            // std::cout<<i <<"\t" << tx_int <<"\t" << ty_int <<"\t" << width <<"\t" << n_rows <<"\t" << n_cols << std::endl;
            //print_transformation_matrix(tr);
            output_image.rowRange(ty_int, ty_int + n_rows).colRange(tx_int, tx_int+n_cols) += tmp;
//            cv::namedWindow( "Display window", cv::WINDOW_AUTOSIZE );// Create a window for display.

//            show_cv_matrix(tmp, "SlitImage");
//            cv::waitKey(5);
            // tmp.copyTo(img.rowRange(ty_int, ty_int + n_rows).colRange(tx_int, tx_int+n_cols));
        }
    }
    return 0;
}

#ifdef USE_GPU
void MatrixSimulator::simulate_spectrum(cv::gpu::GpuMat& slit_image)
{
//    cv::Mat img_cpu = cv::Mat::zeros(4096*3, 4096*3, slit_image.type());
    cv::gpu::GpuMat img = cv::gpu::GpuMat(this->ccd->data);
#pragma omp parallel for
    for(int o=this->orders.front(); o<this->orders.back()+1; ++o){
        std::cout << o << std::endl;
        this->simulate_order(o, slit_image, img, false);
    }
    //img.download(img_cpu);
    //return img_cpu;
}
#endif

void MatrixSimulator::simulate_spectrum(bool aberrations)
{
    this->set_efficiencies(this->efficiencies);

    this->prepare_sources(this->sources);
    //cv::Mat img = cv::Mat::zeros(4096*3, 4096*3, slit_image.type());
    #pragma omp parallel for
    for(int o=this->orders.front(); o<this->orders.back()+1; ++o){
        this->simulate_order(o, this->slit->slit_image, this->ccd->data, aberrations);
    }
    //return img;
}

void MatrixSimulator::set_efficiencies(std::vector<Efficiency *> &efficiencies)
{
    std::cout<<"Set Efficiencies ... " << std::endl;
    for(auto& o : this->orders)
    {
        this->sim_efficiencies.insert(std::pair<int, std::vector<double> > (o, std::vector<double> (this->sim_wavelength[o].size(), 1.0)));

        // std::vector<double> total_eff(this->sim_wavelength[o].size(), 1.0);
        for(auto& e : efficiencies)
        {
            std::vector<double> sim_eff = e->get_efficieny(o, this->sim_wavelength[o]);
            for(int i=0; i<this->sim_efficiencies[o].size(); ++i)
            {
                this->sim_efficiencies[o][i] *= sim_eff[i];
            }
        }
//        std::cout << "Efficiency " << o << std::endl;
        // this->sim_efficiencies[o] = total_eff;
    }
}

void MatrixSimulator::set_order_range(int min_order, int max_order) {
    this->orders.clear();
    for (int i = min_order; i<max_order+1; ++i)
    {
        this->orders.push_back(i);
    }

}

void MatrixSimulator::prepare_sources(std::vector<Source *> sources) {
    std::cout<<"Prepare sources ..." << std::endl;
    this->sim_spectra.clear();
//    #pragma omp parallel for
    for(auto it =this->orders.begin(); it<this->orders.end(); ++it)
    {
        int o = *it;
//        std::cout << "Order " << o << std::endl;
        if ( this->sim_wavelength.count(o) > 0 ){
            this->sim_spectra.insert(std::pair<int, std::vector<double> > (o, std::vector<double>(this->sim_wavelength[o].size())));
            for(auto& s : sources)
            {
                std::vector<double> spectrum = s->get_spectrum(this->sim_wavelength[o]);
                // vectorToFile(spectrum, "../o"+std::to_string(o)+".dat");
                for (int i = 0; i < spectrum.size(); ++i) {
                    this->sim_spectra[o][i] = spectrum[i];
                }
            }
        }
    }

}

void MatrixSimulator::add_efficiency(Efficiency *eff) {
    this->efficiencies.push_back(eff);
}


void MatrixSimulator::set_ccd(CCD * ccd){
    this->ccd = ccd;
}

void MatrixSimulator::set_slit(Slit * slit){
    this->slit = slit;
}

void MatrixSimulator::set_psfs(PSF * psfs){
    this->psfs = psfs;
}

void MatrixSimulator::add_source(Source *src) {
    this->sources.push_back(src);
}

void MatrixSimulator::save_to_hdf(std::string filename, bool downsample, bool bleed, bool overwrite) {
    this->ccd->save_to_hdf(filename, downsample, bleed, overwrite);
}

void MatrixSimulator::save_to_fits(std::string filename, bool downsample, bool bleed, bool overwrite) {
    this->ccd->save_to_fits(filename, downsample, bleed, overwrite);
}

void MatrixSimulator::save_1d_to_fits(std::string filename) {
    CCfits::FITS infile(filename.c_str(), CCfits::Write);
    string hduName("1Dspectrum_Fiber_"+std::to_string(this->fiber_number));
    int n_orders = this->orders.size();

    CCfits::Table* newTable = infile.addTable(hduName,2);
    for(auto const &o: this->orders)
    {
        newTable->addColumn(CCfits::Tdouble, "Order_"+std::to_string(o),this->sim_spectra[o].size(), "");
        newTable->column("Order_"+std::to_string(o)).write(this->sim_spectra[o],1,2);
        newTable->column("Order_"+std::to_string(o)).write(this->sim_wavelength[o],1,1);
//        colName.push_back("Order_"+std::to_string(o));
//        colForm.push_back(std::to_string(this->sim_spectra[o].size())+"D");
//        colUnit.push_back("");
    }

//    CCfits::Table* newTable = infile.addTable().addTable(hduName,n_orders,colName,colForm,colUnit);

}

void MatrixSimulator::transformation_to_file(std::string filename) {
    std::ofstream myfile;
    myfile.open(filename);
    for(auto const &o: this->orders)
    {
        for (auto const & wl: this->sim_wavelength[o]){
//            cv::Mat tm = this->get_transformation_matrix(o, wl);
//            std::vector<double> tmp;
//            tmp = decompose_matrix(tm);
//            myfile << o << ";" << wl << ";" << tmp[0] << ";" << tmp[1] << ";" << tmp[2] << ";" << tmp[3] << ";"<< tmp[4] << ";" << tmp[5] << std::endl;
        }
    }
    myfile.close();
}