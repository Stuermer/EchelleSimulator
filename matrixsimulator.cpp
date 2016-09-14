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
//#include "opencv2/gpu/gpu.hpp"

void print_transformation_matrix(cv::Mat tm){
  std::cout << tm.at<double>(0,0) << "\t" << tm.at<double>(0,1) << "\t" << tm.at<double>(0,2) << "\n";
  std::cout << tm.at<double>(1,0) << "\t" << tm.at<double>(1,1) << "\t" << tm.at<double>(1,2) << std::endl;
  
}

class CSVRow
{
    public:
        std::string const& operator[](std::size_t index) const
        {
            return m_data[index];
        }
        std::size_t size() const
        {
            return m_data.size();
        }
        void readNextRow(std::istream& str)
        {
            std::string         line;
            std::getline(str, line);

            std::stringstream   lineStream(line);
            std::string         cell;

            m_data.clear();
            while(std::getline(lineStream, cell, ';'))
            {
                m_data.push_back(cell);
            }
        }
    private:
        std::vector<std::string>    m_data;
};

std::istream& operator>>(std::istream& str, CSVRow& data)
{
    data.readNextRow(str);
    return str;
}

class CSVIterator
{
    public:
        typedef std::input_iterator_tag     iterator_category;
        typedef CSVRow                      value_type;
        typedef std::size_t                 difference_type;
        typedef CSVRow*                     pointer;
        typedef CSVRow&                     reference;

        CSVIterator(std::istream& str)  :m_str(str.good()?&str:NULL) { ++(*this); }
        CSVIterator()                   :m_str(NULL) {}

        // Pre Increment
        CSVIterator& operator++()               {if (m_str) { if (!((*m_str) >> m_row)){m_str = NULL;}}return *this;}
        // Post increment
        CSVIterator operator++(int)             {CSVIterator    tmp(*this);++(*this);return tmp;}
        CSVRow const& operator*()   const       {return m_row;}
        CSVRow const* operator->()  const       {return &m_row;}

        bool operator==(CSVIterator const& rhs) {return ((this == &rhs) || ((this->m_str == NULL) && (rhs.m_str == NULL)));}
        bool operator!=(CSVIterator const& rhs) {return !((*this) == rhs);}
    private:
        std::istream*       m_str;
        CSVRow              m_row;
};

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
}

/*QwtPlot * MatrixSimulator::plot_transformations(){

    QwtPlot * plot = new QwtPlot();
    QwtPlotCurve * curve = new QwtPlotCurve();

    std::vector<double> x;
    std::vector<double> y;

    for(auto const& a: this->raw_transformations)
    {
        std::cout << a.first;
        if (a.first==88){
            for (auto rt: a.second)
            {
                y.push_back(rt.wavelength);
                x.push_back(rt.decomposed_matrix[0]);
            }
        }
    }

    curve->setSamples(x.data(), y.data(), x.size());
    curve->attach(plot);
    // plot.show();

    return plot;

}*/

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

/*
cv::gpu::GpuMat MatrixSimulator::transform_slit(cv::gpu::GpuMat& slit_image, cv::Mat& transformation_matrix, double weight){
    int n_rows = abs(round(slit_image.rows * 1.5 * transformation_matrix.at<double>(1,1))); 
    int n_cols = abs(round(slit_image.cols * 4. * transformation_matrix.at<double>(0,0)));
    
    cv::Mat tm = transformation_matrix.clone();
    int tx_int = floor(tm.at<double>(0,2));
    int ty_int = floor(tm.at<double>(1,2));
    tm.at<double>(0,2) += -tx_int + n_rows/2.;
    tm.at<double>(1,2) += -ty_int + n_cols * 0.9;
    
    cv::Mat warp_dst_cpu = cv::Mat::zeros( n_rows, n_cols, slit_image.type() );
    cv::gpu::GpuMat warp_dst=cv::gpu::GpuMat();
    warp_dst.upload(warp_dst_cpu);

    cv::gpu::warpAffine(slit_image.clone(), warp_dst, tm, warp_dst.size() );
    return warp_dst;
}
*/

cv::Mat MatrixSimulator::transform_slit(cv::Mat& slit_image, cv::Mat& transformation_matrix, double& weight){
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
    tm.at<double>(0,2) += -tx_int + n_rows/2.;
    tm.at<double>(1,2) += -ty_int + n_cols * 0.9;
    double det = cv::determinant(transformation_matrix.colRange(0,2));
    cv::Mat warp_dst;
    warp_dst = cv::Mat::zeros( n_rows, n_cols, slit_image.type() );
    cv::warpAffine(slit_image.clone()*weight*10./det , warp_dst, tm, warp_dst.size() );
    return warp_dst;
}

/*
int MatrixSimulator::simulate_order(int order, cv::gpu::GpuMat& slit_image, cv::Mat& output_image)
{
  // cv::Mat img = cv::Mat::zeros(4096*3, 4096*3, CV_64FC1);
  int n = this->sim_matrices[order].size();
  for(int i=0; i<this->sim_matrices[order].size(); ++i){
    cv::Mat tr=this->sim_matrices[order][i];
    cv::gpu::GpuMat tmp = cv::gpu::GpuMat();
    cv::Mat tmp2;
    tmp = transform_slit(slit_image, tr, this->sim_efficiencies[order][i]);
    tmp.download(tmp2);
    int tx_int = floor(tr.at<double>(0,2));
    int ty_int = floor(tr.at<double>(1,2));   
    
    int n_rows = abs(round(slit_image.rows * 1.5 * tr.at<double>(1,1))); 
    int n_cols = abs(round(slit_image.cols * 4. * tr.at<double>(0,0)));
    

    //int TSX = -tx_int + n_rows/2.;
    //tm.at<double>(1,2) += -ty_int + n_cols * 0.9;
    
    int left_margin = std::max(0, tx_int);
    int right_margin = std::min(output_image.cols, tx_int+n_cols);
    int width = right_margin - left_margin;

    if (width>0 && width>=n_cols && (ty_int>0) && (ty_int+n_rows<output_image.rows) && tx_int>0 && tx_int+n_cols<output_image.cols)
    {
      // std::cout<<i <<"\t" << tx_int <<"\t" << ty_int <<"\t" << width <<"\t" << n_rows <<"\t" << n_cols << std::endl;
      //print_transformation_matrix(tr);
      output_image.rowRange(ty_int, ty_int + n_rows).colRange(tx_int, tx_int+n_cols) += tmp2;
      // tmp.copyTo(img.rowRange(ty_int, ty_int + n_rows).colRange(tx_int, tx_int+n_cols));  
    }
  }
  return 0;
}
*/

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
          cv::Mat psf = this->psfs->get_PSF(order, wavelength);
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

/*
cv::Mat MatrixSimulator::simulate_spectrum(cv::gpu::GpuMat& slit_image)
{
  cv::Mat img = cv::Mat::zeros(4096*3, 4096*3, slit_image.type());
#pragma omp parallel for
  for(int o=this->orders.front(); o<this->orders.back()+1; ++o){
    std::cout << o << std::endl;
    this->simulate_order(o, slit_image, img);
  }
  return img;
}
*/
cv::Mat MatrixSimulator::simulate_spectrum(cv::Mat& slit_image)
{
  cv::Mat img = cv::Mat::zeros(4096*3, 4096*3, slit_image.type());
#pragma omp parallel for
  for(int o=this->orders.front(); o<this->orders.back()+1; ++o){
    this->simulate_order(o, slit_image, img, true);
  }
  return img;
}

void MatrixSimulator::prepare_efficiencies(std::vector<Efficiency*> &efficiencies)
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

void MatrixSimulator::prepare_sources(std::vector<Source*> sources) {
    std::cout<<"Prepare sources" << std::endl;
    for(auto& o: this->orders)
    {
        std::cout << "Order " << o << std::endl;

        this->sim_spectra.insert(std::pair<int, std::vector<double> > (o, std::vector<double>(this->sim_wavelength[o].size())));
        for(auto& s : sources)
        {
            std::vector<double> spectrum = s->get_spectrum(this->sim_wavelength[o]);
            for(int i=0; i< spectrum.size(); ++i)
            {
                this->sim_spectra[o][i] = spectrum[i];
            }

        }

    }

}


