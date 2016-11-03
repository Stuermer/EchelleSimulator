#include "matrixsimulator.h"
#include <iostream>
#include <fstream>
#include <iterator>
#include "spline.h"
// #include <string>
// #include <stdlib.h>


// using namespace std
#include "helper.h"
#include <opencv2/imgproc.hpp>
#ifdef USE_GPU
#include "opencv2/gpu/gpu.hpp"
#endif

void print_transformation_matrix(cv::Mat tm){
    std::cout << tm.at<double>(0,0) << "\t" << tm.at<double>(0,1) << "\t" << tm.at<double>(0,2) << "\n";
    std::cout << tm.at<double>(1,0) << "\t" << tm.at<double>(1,1) << "\t" << tm.at<double>(1,2) << std::endl;

}



MatrixSimulator::MatrixSimulator()
{
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

        t.transformation_matrix.at<double>(0,0) = std::stod((*loop)[2]);
        t.transformation_matrix.at<double>(0,1) = std::stod((*loop)[3]);
        t.transformation_matrix.at<double>(0,2)= std::stod((*loop)[4]);

        t.transformation_matrix.at<double>(1,0) = std::stod((*loop)[5]);
        t.transformation_matrix.at<double>(1,1) = std::stod((*loop)[6]);
        t.transformation_matrix.at<double>(1,2) = std::stod((*loop)[7]);

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
    std::cout << "Calculate Matrices " << std::endl;
#pragma omp parallel for
    for(int o=this->orders.front(); o<this->orders.back()+1; ++o)
    {
        std::cout<< o << std::endl;
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

cv::Mat MatrixSimulator::get_transformation_matrix(int o, double wavelength){
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

cv::Mat MatrixSimulator::transform_slit(cv::Mat& slit_image, cv::Mat& transformation_matrix, double weight){
    double sx, sy, a,b,c,d;
    a= transformation_matrix.at<double>(0,0);
    b= transformation_matrix.at<double>(1,0);
    c= transformation_matrix.at<double>(0,1);
    d= transformation_matrix.at<double>(1,1);
    sx = sqrt(a*a+b*b);
    sy = sqrt(c*c+d*d);
    int n_rows = abs(round(slit_image.rows * 1.5 * sy));
    int n_cols = abs(round(slit_image.cols * 4. * sx));

    cv::Mat tm = transformation_matrix.clone();
    int tx_int = floor(tm.at<double>(0,2));
    int ty_int = floor(tm.at<double>(1,2));
    tm.at<double>(0,2) += -tx_int + n_cols/2.;
    tm.at<double>(1,2) += -ty_int + n_rows * 0.75;
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

int MatrixSimulator::simulate_order(int order, cv::Mat& slit_image, cv::Mat& output_image, bool aberrations)
{

    // cv::Mat img = cv::Mat::zeros(4096*3, 4096*3, CV_64FC1);
    int n = this->sim_matrices[order].size();
    std::cout << "Simulate order " << order << " number of wavelength " << n << std::endl;

    for(int i=0; i<n; ++i)
    {
        cv::Mat tr = this->sim_matrices[order][i];
        double weight = this->sim_efficiencies[order][i];
        double wavelength = this->sim_wavelength[order][i];
        weight *= this->sim_spectra[order][i];
        cv::Mat tmp = transform_slit(slit_image, tr, weight);

        if (aberrations){
            PSF * psf_p = dynamic_cast<PSF *> (this->psfs);
            cv::Mat psf = psf_p->get_PSF(order, wavelength);
            cv::flip(psf, psf, -1);
            cv::filter2D(tmp, tmp, -1, psf, cvPoint(psf.cols/2, psf.rows/2));
        }

        int tx_int = floor(tr.at<double>(0,2));
        int ty_int = floor(tr.at<double>(1,2));

        double sx, sy, a,b,c,d;
        a= tr.at<double>(0,0);
        b= tr.at<double>(1,0);
        c= tr.at<double>(0,1);
        d= tr.at<double>(1,1);
        sx = sqrt(a*a+b*b);
        sy = sqrt(c*c+d*d);

        int n_rows = abs(round(slit_image.rows * 1.5 * sy));
        int n_cols = abs(round(slit_image.cols * 4. * sx));


        //int TSX = -tx_int + n_rows/2.;
        //tm.at<double>(1,2) += -ty_int + n_cols * 0.9;

        int left_margin = std::max(0, tx_int);
        int right_margin = std::min(output_image.cols, tx_int+n_cols);
        int width = right_margin - left_margin;

        if (width>0 && width>=n_cols && (ty_int>0) && (ty_int+n_rows<output_image.rows) && tx_int>0 && tx_int+n_cols<output_image.cols)
        {
            // std::cout<<i <<"\t" << tx_int <<"\t" << ty_int <<"\t" << width <<"\t" << n_rows <<"\t" << n_cols << std::endl;
            //print_transformation_matrix(tr);
            output_image.rowRange(ty_int, ty_int + n_rows).colRange(tx_int, tx_int+n_cols) += tmp;
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

void MatrixSimulator::simulate_spectrum()
{
    this->set_efficiencies(this->efficiencies);

    this->prepare_sources(this->sources);
    //cv::Mat img = cv::Mat::zeros(4096*3, 4096*3, slit_image.type());
#pragma omp parallel for
    for(int o=this->orders.front(); o<this->orders.back()+1; ++o){
        this->simulate_order(o, this->slit->slit_image, this->ccd->data, true);
    }
    //return img;
}

void MatrixSimulator::set_efficiencies(std::vector<Efficiency *> &efficiencies)
{
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
        std::cout << "Efficiency " << o << std::endl;
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
    std::cout<<"Prepare sources" << std::endl;
    this->sim_spectra.clear();
    for(auto const& o: this->orders)
    {
        std::cout << "Order " << o << std::endl;
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

void MatrixSimulator::save_to_file(std::string filename, bool downsample, bool bleed, bool overwrite) {
    this->ccd->save_to_file(filename, downsample, bleed, overwrite);
}

void MatrixSimulator::transformation_to_file(std::string filename) {
    std::ofstream myfile;
    myfile.open(filename);
    for(auto const &o: this->orders)
    {
        for (auto const & wl: this->sim_wavelength[o]){
            cv::Mat tm = this->get_transformation_matrix(o, wl);
            std::vector<double> tmp;
            tmp = decompose_matrix(tm);
            myfile << o << ";" << wl << ";" << tmp[0] << ";" << tmp[1] << ";" << tmp[2] << ";" << tmp[3] << ";"<< tmp[4] << ";" << tmp[5] << std::endl;
        }
    }
    myfile.close();
}