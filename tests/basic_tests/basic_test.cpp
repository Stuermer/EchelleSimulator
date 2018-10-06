#include "gtest/gtest.h"
#include "helper.h"


TEST(basic_test, test_phoenix) {
    EXPECT_EQ(0, download_wave_grid("tmp.fits"));
    EXPECT_EQ(1, check_for_file("tmp.fits"));
    remove("tmp.fits");
    EXPECT_EQ(0, check_for_file("tmp.fits"));

    // try invalid parameters
    EXPECT_EQ(1, download_phoenix(0, 0, 0, 0, "tmp.fits"));
    // try valid parameters
    EXPECT_EQ(0, download_phoenix(3200, 5.5, 0., 0., "tmp.fits"));
    EXPECT_EQ(1, check_for_file("tmp.fits"));
    remove("tmp.fits");
    EXPECT_EQ(0, check_for_file("tmp.fits"));
}

TEST(basic_test, test_interpolation) {
    std::cout<<"Test interpolation" << std::endl;
    std::map<double, double> x = { {0., 0.}, {1., 1.}};
    EXPECT_EQ(0.5, interpolate(x, 0.5));
}

// TODO: Write tests for tracing, sources, random_generator etc.