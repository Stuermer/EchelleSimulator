#include "CCD.h"
#include <CCfits/CCfits>
#include "helper.h"

CCD::CCD(int Nx, int Ny, double pixelsize) : Nx(Nx), Ny(Ny), pixel_size(pixelsize) {
    {
        this->data = std::vector<int>(Ny * Nx, 0);
    }
}

double *CCD::get_pixelsize() {
    return &this->pixel_size;
}

void CCD::save_to_hdf(std::string filename, bool down_sample, bool bleed, bool overwrite) {
//    cv::Mat res = this->get_image(down_sample, bleed);
//    hdf5opencv::hdf5save(filename.c_str(), "image", res, overwrite);
    // MatToFile(res, filename);
}

void CCD::save_to_fits(std::string filename, bool overwrite) {
    int n_axis = 2;
    long n_axes[2] = {Nx, Ny};
    std::unique_ptr<CCfits::FITS> pFits;
    // Try to read in old image
    long old_img_ax1 = 0;
    long old_img_ax2 = 0;
    std::valarray<int> contents;
    try {

        pFits.reset(new CCfits::FITS(filename, CCfits::Read));
        CCfits::PHDU &image = pFits->pHDU();

        image.readAllKeys();
        image.read(contents);

        old_img_ax1 = image.axis(0);
        old_img_ax2 = image.axis(1);
        pFits->destroy();
        std::cout << "Old data found" << std::endl;
    }
    catch (...) {
        // no old data
        std::cout << "No old data found" << std::endl;;
        pFits.reset(new CCfits::FITS(filename, LONG_IMG, n_axis, n_axes));
        pFits->flush();
        pFits->destroy();
    }

    try {
        pFits.reset(new CCfits::FITS(filename, CCfits::Write));

    }
    catch (CCfits::FITS::CantCreate) {
        // ... or not, as the case may be.
        std::cout << "Can't create FITS file." << std::endl;
    }
    long n_elements(1);
    n_elements = std::accumulate(&n_axes[0], &n_axes[n_axis], 1, std::multiplies<long>());
    long f_pixel(1);
    std::valarray<int> data_array(Nx * Ny);
    if ((old_img_ax1 > 0) & (old_img_ax2 > 0) & !overwrite) {
        for (int i = 0; i < Nx * Ny; ++i) {
            data_array[i] = (int) this->data[i] + contents[i];
        }
    } else {

        for (int i = 0; i < Nx * Ny; ++i) {
            data_array[i] = (int) this->data[i];
        }
    }

    pFits->pHDU().write(f_pixel, n_elements, data_array);
    pFits->flush();

}

CCD::~CCD() = default;

//cv::Mat CCD::get_image(bool down_sample, bool bleed) {
//    cv::Mat result;
//    cv::Mat initial(Ny,Nx, CV_16UC1);
//    for(int i=0; i<Ny; ++i){
//        for (int j = 0; j < Nx; ++j) {
//            initial.at<unsigned short>(j, i) = (unsigned short) this->data[j*Nx+i];
//        }
//    }
//    if (down_sample){
//        #ifdef  USE_GPU
//                {
//                this->data.download(result)
//            };
//        #else
//                {
//                    result = cv::Mat(Ny, Nx, initial.type());
//                }
//        #endif
//        cv::resize(initial, result, result.size(), cv::INTER_NEAREST);
//    }
//    else {
//        result = initial;
//    }
//
//    if (bleed){
//        this->do_bleed(result, 650000.);
//    }
//
//    return result;
//}

//void do_bleed_up_down(cv::Mat &input, double & limit, int i, int j){
//    //find limit
//    int k = i;
//    while ((k>1) && (k<input.rows-1) && (input.at<unsigned short>(k, j)))
//    {
//        k++;
//    }
//
////    cv::blur(input.colRange(i-1,k+1), input.colRange(i-1,k+1), cv::Size(1,3),cv::Point(0.,1.5) ,cv::BORDER_REPLICATE);
//    if(i>0)
//    {
//        if (input.at<unsigned short>(i-1,j) > limit)
//            do_bleed_up_down(input, limit, i-1, j);
//    }
//    double diff = input.at<double>(i, j) - limit;
//    if ((diff > 0) && (i + 2 < input.rows) && (i - 2 > 0)) {
//        std::cout<< i << "\t" << j << std::endl;
//        input.at<double>(i + 1, j) += diff/2.;
//        input.at<double>(i - 1, j) += diff/2.;
////        input.at<double>(i - 1, j) += diff / 2.;
//        input.at<double>(i, j) -= diff;
//        do_bleed_up_down(input, limit, i - 1, j);
//        do_bleed_up_down(input, limit, i + 1, j);
//    }
//}

//void do_bleed_down(cv::Mat &input, double limit, int i, int j){
//    double diff = input.at<double>(i,j) - limit;
//    if ((diff>0) && (i+1 < input.rows) && (i-1 >0))  {
//        input.at<double>(i - 1, j) += diff / 2.;
//        do_bleed_down(input, limit, i - 1, j);
//    }
//
//}

//void CCD::do_bleed(cv::Mat &input, double limit) {
//
//    for(int i=0; i<input.cols; ++i) {
//        for (int j = 0; j < input.rows; ++j) {
//            if(input.at<unsigned short>(j,i) > limit){
//                do_bleed_up_down(input, limit, j, i);
//            }
//        }
//    }
//}
//
