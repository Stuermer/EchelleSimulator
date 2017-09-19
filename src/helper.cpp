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
#include <curl/curl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string>

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

//void show_cv_matrix(cv::Mat img, std::string windowname="image") {
//    double minVal, maxVal;
//
//    cv::Mat img_show = img.clone();
//    cv::minMaxLoc(img_show, &minVal, &maxVal); //find minimum and maximum intensities
//    int ty = img_show.type();
//    img_show.convertTo(img_show,CV_8U,255.0/(maxVal - minVal), -minVal * 255.0/(maxVal - minVal));
//
//     cv::cvtColor(img_show, img_show, CV_GRAY2RGB);
//
//    cv::namedWindow(windowname,CV_WINDOW_NORMAL);
//    cv::imshow(windowname, img_show);
//    cv::resizeWindow(windowname, 1024,1024);
//    cv::waitKey(1);
//
//}

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

herr_t file_info(hid_t loc_id, const char *name, const H5L_info_t *linfo, void *opdata)
{
    hid_t group;
    auto group_names=reinterpret_cast< std::vector<std::string>* >(opdata);
    group_names->push_back(name);
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

size_t write_data(void *ptr, size_t size, size_t nmemb, FILE *stream) {
    size_t written = fwrite(ptr, size, nmemb, stream);
    return written;
}

int download_phoenix(std::string Teff, std::string log_g, std::string z, std::string alpha){

    //input positive z without '+'
    if(z.at(0) != '-' && z.at(0) != '+' && z.at(0) != '0'){

        z = "+" + z;

    }

    //input positive alpha without '+'
    if(alpha.at(0) != '-' && alpha.at(0) != '+' && alpha.at(0) != '0'){

        alpha = "+" + alpha;

    }

    //fix attempts at inputting z = -0.0
    if(z == "0"){

        z = "-" + z + ".0";

    }
    else if(z == "0."){

        z = "-" + z + "0";

    }
    else if(z == "0.0"){

        z = "-" + z;

    }
    else{

    }

    //fix attempts at inputting alpha = -0.0
    if(alpha == "0"){

        alpha = "-" + alpha + ".00";

    }
    else if(alpha == "0."){

        alpha = "-" + alpha + "00";

    }
    else if(alpha == "0.0"){

        alpha = "-" + alpha + "0";

    }
    else if(alpha == "0.00"){

        alpha = "-" + alpha;

    }
    else{

    }

    //fix attempts at inputting log_g = -0.00
    if(log_g == "0"){

        log_g = "-" + log_g + ".00";

    }
    else if(log_g == "0."){

        log_g = "-" + log_g + "00";

    }
    else if(log_g == "0.0"){

        log_g = "-" + log_g + "0";

    }
    else if(log_g == "0.00"){

        log_g = "-" + log_g;

    }
    else{

    }

    //pad for z inputs
    if(z.length() == 2){
        z = z + ".0";
    }
    else if(z.length() == 3){
        z = z + "0";
    }
    else{

    }

    //pad for log_g inputs
    if(log_g.length() == 1){
        log_g = log_g + ".00";
    }
    else if(log_g.length() == 2){
        log_g = log_g + "00";
    }
    else if(log_g.length() == 3){
        log_g = log_g + "0";
    }
    else{

    }

    //pad for alpha inputs
    if(alpha.length() == 2){
        alpha = alpha + ".00";
    }
    else if(alpha.length() == 3){
        alpha = alpha + "00";
    }
    else if(alpha.length() == 4){
        alpha = alpha + "0";
    }

    std::string url = "ftp://phoenix.astro.physik.uni-goettingen.de/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/";

    url += "Z" + z;

    if(alpha != "-0.00"){

        url += ".Alpha=" + alpha;

    }

    url += "/lte";

    for(int i = 0; i < 5 - Teff.length(); i++){

        Teff = "0" + Teff;

    }

    url += Teff + "-" + log_g;

    url += z;

    if(alpha != "-0.00"){

        url += ".Alpha=" + alpha;

    }

    url += ".PHOENIX-ACES-AGSS-COND-2011-HiRes.fits";

    std::cout<<"Downloading spectra from: "<<url<<std::endl; //cover z = 0 case also see if adding x.0 and + | - is doable

    CURL *curl;
    FILE *fp;
    CURLcode res;
    std::string path = "../data/phoenix_spectra/test.fits";
    remove(path.c_str());

    curl = curl_easy_init();
    if (curl) {
        fp = fopen(path.c_str(), "wb");
        curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_data);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, fp);
        res = curl_easy_perform(curl);
        /* always cleanup */
        curl_easy_cleanup(curl);
        fclose(fp);
    }

    return res;


}

int download_wave_grid(std::string path){

    std::string url = "ftp://phoenix.astro.physik.uni-goettingen.de/HiResFITS//WAVE_PHOENIX-ACES-AGSS-COND-2011.fits";

    CURL *curl;
    FILE *fp;
    CURLcode res;

    curl = curl_easy_init();
    if (curl) {
        fp = fopen(path.c_str(), "wb");
        curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_data);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, fp);
        res = curl_easy_perform(curl);
        /* always cleanup */
        curl_easy_cleanup(curl);
        fclose(fp);
    }

    return res;

}

bool check_for (const std::string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }
}