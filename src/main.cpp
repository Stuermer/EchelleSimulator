#include <opencv2/imgproc.hpp>
#include <iostream>

#include "matrixsimulator.h"
#include "slit.h"
#include "efficiency.h"
#include <chrono>
#include "helper.h"
#include "noise.h"
#include "source.h"
#include "PSF.h"

using namespace std::chrono;

int main(int argc, char *argv[])
{
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    PhoenixSpectrum ps = PhoenixSpectrum("/home/julian/Dissertation/CRIRES-POP/template/lte03400-4.00-1.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits","/home/julian/Dissertation/CRIRES-POP/template/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits", 0.45, 0.7);

    PSF_ZEMAX psfs = PSF_ZEMAX(argv[1]);
//    PSF_gaussian psfs = PSF_gaussian(1., 10.);
//    cv::Mat r = psfs.get_PSF(10,3.);
//    print_cv_matrix_info(r, "psf info");
//
    MatrixSimulator m;
    m.psfs = &psfs;
    m.read_transformations(argv[2]);
    //m.set_order_range(89,95);
    m.calc_splines();
    m.set_wavelength(10000);
    m.calc_sim_matrices();
    slit s = slit(50., 150., 10);
    std::cout <<" calc spectrum" << std::endl;
    std::vector<double> spec = ps.get_spectrum(m.sim_wavelength[89]);
    std::cout << "done " << std::endl;
    vectorToFile(spec, "../1dspec");

    GratingEfficiency ge = GratingEfficiency(0.8, 76., 76., 31.6);
    std::vector<Efficiency*> efficiencies;
    efficiencies.push_back(&ge);
    m.prepare_efficiencies(efficiencies);

    IdealEtalon cs = IdealEtalon(10., 1., 0., 0.9);
    std::vector<Source*> sources;
    sources.push_back(&cs);
    m.prepare_sources(sources);

    cv::Mat img = m.simulate_spectrum(s.slit_image);

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>( t2 - t1 ).count()/1000000.;
    std::cout << "Duration: "  << duration << std::endl;

    cv::Mat img_ccd;
    img_ccd = cv::Mat::zeros(4096, 4096, s.slit_image.type());

    cv::resize(img, img_ccd, img_ccd.size(), cv::INTER_NEAREST);


    double minVal, maxVal;
    cv::minMaxLoc(img_ccd, &minVal, &maxVal); //find minimum and maximum intensities

    std::cout << minVal << " min and max " << maxVal <<std::endl;

    //noise n = noise();
    //n.add_poisson(img_ccd);

    MatToFile(img_ccd, "../image2.dat");

    //img_ccd.convertTo(img_ccd, CV_8UC1, 255/(maxVal - minVal), -minVal * 255/(maxVal - minVal));
    //cv::minMaxLoc(img_ccd, &minVal, &maxVal); //find minimum and maximum intensities
    //std::cout << minVal << " min and max " << maxVal <<std::endl;


    // cv::cvtColor(img_ccd, img_ccd, CV_GRAY2RGB);
    // cv::imwrite("../test.png", img_ccd );

   // cv::namedWindow("image", CV_WINDOW_AUTOSIZE);

    //cv::imshow("image", img_ccd);
    //cv::waitKey(0);
   /// Wait until user exits the program
  //  topLevelLabel.show();
  /*
    QLabel topLevelLabel;
    QPixmap pixmap = QPixmap::fromImage(QtOcv::mat2Image(img_ccd));

    topLevelLabel.setPixmap(pixmap);
    topLevelLabel.setMask(pixmap.mask());
    topLevelLabel.show();*/
    // return app.exec();
    return 0;
}
