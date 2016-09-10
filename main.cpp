#include <QtWidgets>
#include <QtGui>
#include <opencv2/imgproc.hpp>
#include <iostream>

#include "matrixsimulator.h"
#include "slit.h"
#include "efficiency.h"
#include "cvmatandqimage.h"
#include <chrono>

using namespace std::chrono;

int main(int argc, char *argv[])
{
//    MainWindow w;
//    w.show();
    
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    QApplication app(argc, argv);

    MatrixSimulator m;
    m.read_transformations("/home/julian/Repos/Python_projects/PyEchelle-github/PyEchelleSimulator/test.csv");
    m.calc_splines();
    m.set_wavelength(10000);
    m.calc_sim_matrices();
    slit s = slit(50., 150., 10);

    GratingEfficiency ge = GratingEfficiency();
    std::vector<GratingEfficiency> efficiencies;
    efficiencies.push_back(ge);

    
    m.prepare_efficienies(efficiencies);
    
    //s.show();
    // s.show();
    cv::Mat img = m.simulate_spectrum(s.slit_image);
    // cv::Mat img = cv::Mat::zeros(4096*3, 4096*3, s.slit_image.type());
    // int i = m.simulate_order(119, s.slit_image, img);
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>( t2 - t1 ).count()/1000000.;
    std::cout << "Duration: "  << duration << std::endl;
    cv::Mat img_ccd = cv::Mat::zeros(1024, 1024, s.slit_image.type());
    cv::resize(img, img_ccd, img_ccd.size());
    cv::Mat draw;
    /*
    double minVal, maxVal;
    cv::minMaxLoc(img_ccd, &minVal, &maxVal); //find minimum and maximum intensities
    img_ccd /= maxVal;
//    std::cout << minVal << " min and max " << maxVal <<std::endl;
//    img_ccd.convertTo(draw, CV_8UC1, 255/(maxVal - minVal), -minVal * 255/(maxVal - minVal));
    cv::minMaxLoc(img_ccd, &minVal, &maxVal); //find minimum and maximum intensities
    std::cout << minVal << " min and max " << maxVal <<std::endl;
    cv::namedWindow("image", CV_WINDOW_AUTOSIZE);
    cv::imshow("image", img_ccd);
    cv::waitKey(0);
   /// Wait until user exits the program
  //  topLevelLabel.show();
  
    QLabel topLevelLabel;
    QPixmap pixmap = QPixmap::fromImage(QtOcv::mat2Image(img_ccd));
        
    topLevelLabel.setPixmap(pixmap);
    topLevelLabel.setMask(pixmap.mask());
    topLevelLabel.show();*/
    return app.exec();
    // return 0;
}
