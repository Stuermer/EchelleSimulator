#ifndef SLIT_H
#define SLIT_H
#include <opencv2/core.hpp>
#ifdef USE_GPU
    #include "opencv2/gpu/gpu.hpp"
#endif

/**
 * \class Slit
 * \brief Class handles the spectrograph input slit
 *
 * This class represents the spectrograph intput slit.
 * The slit has a certain width and height (usually given by the spectrograph model) and is represented by a 2D array.
 *
 * The intensity distribution within the slit can be used to mimic different slit shapes e.g. a circular fiber,
 * a rectangular fiber, an octagonal fiber by filling the intensity distribution of the slit accordingly.
 * It can also be used to simulate a sliced fiber.
 *
 * As the input slit is represented by an image, it could in principle also be used to simulate intensity variations
 * within the slit to mimic incomplete fiber scrambling.
 */
class Slit
{
public:
    /**
     * Default constructor
     */
    Slit();

    /**
     * Constructor for use with spectrograph model
     * @param w width in microns
     * @param h height in microns
     * @param slit_sampling number of pixels used for sampling (given by spectrograph model)
     */
    Slit(double w, double h, int slit_sampling);

    /**
     * Sets a slit to a certain width, height and sampling. Usually given by spectrograph model.
     * @param w width in microns
     * @param h heigh in microns
     * @param slit_sampling number of pixels used for sampling the slit.
     */
    void set_slit(double w, double h, int slit_sampling);

    /**
     * Plots the slit.
     */
    void show();

//private:
    double w;
    double h;
    int w_px;
    int h_px;
    double ratio;

    int slit_sampling;

    #ifdef USE_GPU
        cv::gpu::GpuMat slit_image;
        bool use_gpu = true;
    #else
        cv::Mat slit_image;
        bool use_gpu = false;
    #endif
};

#endif // SLIT_H
