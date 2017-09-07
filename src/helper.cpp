#include "helper.h"

#include <vector>
#include <cmath>
#include "opencv2/core.hpp"

#include <vector>
#include <iterator>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <opencv2/imgproc.hpp>
#include <highgui.h>
//#include <CCfits>
#include <CCfits/FITS.h>
#include <CCfits/ExtHDU.h>
#include <Eigen/Dense>
#include <map>

void vectorToFile(std::vector<double> const& vec, std::string const& filename) {
  std::ofstream file(filename);
  auto first = true;
  for (float f : vec) 
  { 
      if (!first) { file << ","; } 
      first = false; 
      file << f; 
  }
  file << std::endl;
  file.close();  
}

std::vector<double> decompose_matrix(Matrix23f mat){
    std::vector<double> result;
    /*
     * Matrix looks like:
     * a b tx
     * d e ty
     */

    double a = mat(0,0);
    double b = mat(0,1);
    double d = mat(1,0);
    double e = mat(1,1);
    double tx = mat(0,2);
    double ty = mat(1,2);
    
    double sx = sqrt(a*a+d*d);
    double sy = sqrt(b*b+e*e);
    
    double phi = atan2(d,a);
    if (phi<0.1)
      phi += 2.*M_PI;
    
    double shear = atan2(-b,e) - phi;
    if (shear < -6.1)
        shear += 2.*M_PI;

    result.push_back(sx);
    result.push_back(sy);
    result.push_back(shear);
    result.push_back(phi);
    result.push_back(tx);
    result.push_back(ty);
    // return <sx, sy, shear, rot, tx ,ty>
    return result;
}

Matrix23f compose_matrix(std::vector<double> parameters){
  double sx = parameters[0];
  double sy = parameters[1];
  double shear = parameters[2];
  double rot = parameters[3];

  Matrix23f m;
    m(0,0) = sx*cos(rot);
    m(1,0) = sx*sin(rot);
    m(0,1) = -sy*sin(rot+shear);
    m(1,1)= sy*cos(rot+shear);
    m(0,2) = parameters[4];
    m(1,2) = parameters[5];
  return m;
}

std::vector<std::size_t> compute_sort_order(const std::vector<double> &v) {
    std::vector<std::size_t> indices(v.size());
    std::iota(indices.begin(), indices.end(), 0u);
    std::sort(indices.begin(), indices.end(), [&](int lhs, int rhs) {
        return v[lhs] < v[rhs];
    });
    std::vector<std::size_t> res(v.size());
    for (std::size_t i = 0; i != indices.size(); ++i) {
        res[indices[i]] = i;
    }
    return res;
}

void show_cv_matrix(cv::Mat img, std::string windowname="image") {
//    double minVal, maxVal;

  //  cv::Mat img_show = img.clone();

//    cv::Mat img_show = cv::Mat::zeros(512/3, 1024, img.type());
//    cv::resize(img.rowRange(0, 512*3).colRange(0,4096*3), img_show, img_show.size(), cv::INTER_NEAREST);
//
//    cv::minMaxLoc(img_show, &minVal, &maxVal); //find minimum and maximum intensities
//    int ty = img_show.type();
    //img_show.convertTo(img_show,CV_8U,255.0/(maxVal - minVal), -minVal * 255.0/(maxVal - minVal));

    // cv::cvtColor(img_show, img_show, CV_GRAY2RGB);

//    cv::namedWindow(windowname,CV_WINDOW_NORMAL);
  //  cv::imshow(windowname, img_show);
   // cv::resizeWindow(windowname, 1024,1024);
   // cv::waitKey(1);

}

void print_cv_matrix_info(cv::Mat img, std::string imagename="Image") {
    std::cout<< "------------" << imagename << "----------------:" << std::endl;
    std::cout << "Dimensions: " << img.rows << " x " << img.cols << std::endl;
    std::cout << "Type: " << img.type() << std::endl;
    double minVal, maxVal;
    cv::minMaxLoc(img, &minVal, &maxVal); //find minimum and maximum intensities
    std::cout << "Min/Max: " << minVal << " / " << maxVal << std::endl;
    std::cout<< "------------------------------------------:"<<std::endl << std::endl;


}

double wrap_rads(double r)
{
    while ( r > M_PI ) {
        r -= 2 * M_PI;
    }

    while ( r <= -M_PI ) {
        r += 2 * M_PI;
    }

    return r;
}

void create_fits_file(std::string filename){

    CCfits::FITS infile(filename.c_str(), CCfits::Write);

};

double interpolate(const std::map<double,double> &data,
                    double x)
{
    typedef std::map<double, double>::const_iterator i_t;

    i_t i=data.upper_bound(x);
    if(i==data.end())
    {
        return (--i)->second;
    }
    if (i==data.begin())
    {
        return i->second;
    }
    i_t l=i; --l;

    const double delta=(x- l->first)/(i->first - l->first);
    return delta*i->second +(1-delta)*l->second;
}

//int save_to_hdf(const std::string filename, cv::Mat img){
//    std::auto_ptr<CCfits::FITS> pFits(0);
//
//    long naxis    =   2;
//    long naxes[2] = { img.cols, img.rows };
//
//    try
//    {
//        pFits.reset( new CCfits::FITS(filename , DOUBLE_IMG , naxis , naxes ) );
//    }
//    catch (CCfits::FITS::CantOpen)
//    {
//        return -1;
//    }
//
//    long& vectorLength = naxes[0];
//    long& numberOfRows = naxes[1];
//    long nelements(1);
//
//
//    // Find the total size of the array.
//    // this is a little fancier than necessary ( It's only
//    // calculating naxes[0]*naxes[1]) but it demonstrates  use of the
//    // C++ standard library accumulate algorithm.
//
//    nelements = std::accumulate(&naxes[0],&naxes[naxis],1,std::multiplies<long>());
//
//    std::vector<long> extAx ;
//    extAx.push_back(img.rows);
//    extAx.push_back(img.cols);
//
//    string newName ("IMAGE");
//    CCfits::ExtHDU* imageExt = pFits->addImage(newName,FLOAT_IMG, extAx);
//
//
//    std::valarray<double> array(nelements);
//    for (int i = 0; i < numberOfRows; ++i)
//    {
//        for (int j =0; j<img.cols; ++j)
//            array[i*numberOfRows+j] = img.at<double>(i, j);
//    }
//
//    long  fpixel(1);
//
//    imageExt->write(fpixel, (long) img.cols, array);
////    imageExt->write(fpixel, img.cols, array);
//
//    return 0;
//
//};

herr_t file_info(hid_t loc_id, const char *name, const H5L_info_t *linfo, void *opdata)
{
    hid_t group;
    auto group_names=reinterpret_cast< std::vector<std::string>* >(opdata);
//    group = H5Gopen2(loc_id, name, H5P_DEFAULT);

    group_names->push_back(name);
//    std::cout << "Name : " << name << std::endl;
//    H5Gclose(group);
    return 0;
}

std::vector<float> random_from_2_distributions(std::vector<float> wl, std::vector<float> density1, std::vector<float> density2, int N_samples){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::vector<float> result;
    std::vector<float> combined_vec;
    for(int i=0; i<density1.size(); ++i)
        density1[i] *= density2[i];

    std::piecewise_linear_distribution<> combined_dis(wl.begin(), wl.end(), combined_vec.begin());
    for(int i=0; i<N_samples; ++i)
        result.push_back(combined_dis(gen));
    return result;
}

