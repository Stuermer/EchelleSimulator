#include "helper.h"
#include <vector>
#include <cmath>
#include <vector>
#include <iterator>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <CCfits/FITS.h>
#include <CCfits/ExtHDU.h>
#include <map>
#include <curl/curl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string>
/// mark fmt as header only
#define FMT_HEADER_ONLY
#include <fmt/format.h>

void vector_to_file(std::vector<double> const &vec, std::string const &filename) {
    std::ofstream file(filename);
    auto first = true;
    for (float f : vec) {
        if (!first) { file << ","; }
        first = false;
        file << f;
    }
    file << std::endl;
    file.close();
}

std::array<double, 6> decompose_matrix(std::array<double, 6> mat) {

    /*
     * Matrix looks like:
     * a b tx   =   m0 m1 m2
     * d e ty       m3 m4 m5
     */

    double a = mat[0];
    double b = mat[1];
    double d = mat[3];
    double e = mat[4];
    double tx = mat[2];
    double ty = mat[5];

    double sx = sqrt(a * a + d * d);
    double sy = sqrt(b * b + e * e);

    double phi = atan2(d, a);
    if (phi < 0.1)
        phi += 2. * M_PI;

    double shear = atan2(-b, e) - phi;
    if (shear < -6.1)
        shear += 2. * M_PI;

    std::array<double, 6> result = {sx, sy, shear, phi, tx, ty};
    // return <sx, sy, shear, rot, tx ,ty>
    return result;
}

std::array<double, 6> compose_matrix(std::vector<double> parameters) {
    double sx = parameters[0];
    double sy = parameters[1];
    double shear = parameters[2];
    double rot = parameters[3];

    std::array<double, 6> m = {
            sx * cos(rot),
            -sy * sin(rot + shear),
            parameters[4],
            sx * sin(rot),
            sy * cos(rot + shear),
            parameters[5],
    };
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


double wrap_rads(double r) {
    while (r > M_PI) {
        r -= 2 * M_PI;
    }

    while (r <= -M_PI) {
        r += 2 * M_PI;
    }

    return r;
}

double interpolate(const std::map<double, double> &data,
                   double x) {
    typedef std::map<double, double>::const_iterator i_t;

    i_t i = data.upper_bound(x);
    if (i == data.end()) {
        return (--i)->second;
    }
    if (i == data.begin()) {
        return i->second;
    }
    i_t l = i;
    --l;

    const double delta = (x - l->first) / (i->first - l->first);
    return delta * i->second + (1 - delta) * l->second;
}

herr_t file_info(hid_t loc_id, const char *name, const H5L_info_t *linfo, void *opdata) {
    auto group_names = reinterpret_cast< std::vector<std::string> * >(opdata);
    group_names->push_back(name);
    return 0;
}

size_t write_data(void *ptr, size_t size, size_t nmemb, FILE *stream) {
    size_t written = fwrite(ptr, size, nmemb, stream);
    return written;
}

int download_phoenix(int t_eff, double log_g, double z, double alpha, const std::string path) {
    // check for valid inputs:
    std::cout << "Trying to download PHOENIX Spectrum with: T_eff=" << t_eff << " log_g=" << log_g << " z=" << z
              << " Alpha=" << alpha << std::endl;

    std::vector<int> valid_T;
    std::vector<double> valid_g, valid_z, valid_a;
    for (int i = 2300; i < 7000; i += 100) { valid_T.push_back(i); }
    for (int i = 7000; i < 12200; i += 200) { valid_T.push_back(i); }

    for (int i = 0; i < 12; ++i) { valid_g.push_back(i * 0.5); }
    for (int i = -4; i < -2; ++i) { valid_z.push_back(i); }
    for (int i = -4; i < 3; ++i) { valid_z.push_back(i * 0.5); }

    for (int i = -1; i < 7; ++i) { valid_a.push_back(i * 0.5); }

    bool found = (std::find(valid_T.begin(), valid_T.end(), t_eff) != valid_T.end());
    if (!found) {
        std::cout << "Invalid Effective temperature value for Phoenix spectra." << std::endl;
        return true;
    }

    found = (std::find(valid_g.begin(), valid_g.end(), log_g) != valid_g.end());
    if (!found) {
        std::cout << "Invalid log_g value for Phoenix spectra." << std::endl;
        return true;
    }

    found = (std::find(valid_z.begin(), valid_z.end(), z) != valid_z.end());
    if (!found) {
        std::cout << "Invalid z value for Phoenix spectra." << std::endl;
        return true;
    }

    found = (std::find(valid_a.begin(), valid_a.end(), alpha) != valid_a.end());
    if (!found) {
        std::cout << "Invalid alpha value for Phoenix spectra." << std::endl;
        return true;
    }

    std::string baseurl = "ftp://phoenix.astro.physik.uni-goettingen.de/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/";
    std::string subgrid = "Z";
    (z > 0) ? subgrid += fmt::format("{:+2.1f}", z) : subgrid += fmt::format("{:-2.1f}", z);
    (fabs(alpha) < 1E-10) ? subgrid += "" : subgrid += ".Alpha=";
    (fabs(alpha) < 1E-10) ? subgrid += "" : subgrid += fmt::format("{:+2.2f}", alpha);

    std::string url =
            baseurl + subgrid + "/lte" + fmt::format("{:05}", t_eff) + "-" + fmt::format("{:2.2f}", log_g);
    (z > 0) ? url += fmt::format("{:+2.1f}", z) : url += fmt::format("{:-2.1f}", z);
    (fabs(alpha) < 1E-10) ? url += "" : url += ".Alpha=";
    (fabs(alpha) < 1E-10) ? url += "" : url += fmt::format("{:+2.2f}", alpha);
    url += ".PHOENIX-ACES-AGSS-COND-2011-HiRes.fits";

    std::cout << "Downloading spectra from: " << url
              << std::endl; //cover z = 0 case also see if adding x.0 and + | - is doable

    CURL *curl;
    FILE *fp;
    CURLcode res;
    res=CURLE_OK;

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

int download_wave_grid(std::string path) {

    std::string url = "ftp://phoenix.astro.physik.uni-goettingen.de/HiResFITS//WAVE_PHOENIX-ACES-AGSS-COND-2011.fits";

    CURL *curl;
    FILE *fp;
    CURLcode res;
    res=CURLE_OK;

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

bool check_for_file(const std::string &path) {
    if (FILE *file = fopen(path.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }
}


template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}

std::vector<std::string> split_to_vector(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}