//
// Created by stuermer on 3/27/17.
//

#ifndef ECHELLESIMULATOR_RANDOM_GENERATOR_H
#define ECHELLESIMULATOR_RANDOM_GENERATOR_H
#include <vector>
#include "stdlib.h"
#include <random>
#include <algorithm>

#ifdef USE_CUDA
#include <curand.h>
// Utilities and system includes
#include <helper_functions.h>  // helper for shared functions common to CUDA Samples
#include <helper_cuda.h>       // helper for CUDA Error handling

#endif

template <class T> class RG_uniform_real{
public:
    RG_uniform_real()
    {
        this->distribution = std::uniform_real_distribution<T>();
        std::random_device rdev{};
        this->generator.seed(rdev());
    };


    RG_uniform_real(int seed){
        this->distribution = std::uniform_real_distribution<T>();
        this->generator.seed(seed);
    };

    RG_uniform_real(T min, T max)
    {
        this->distribution = std::uniform_real_distribution<T>(min, max);
        std::random_device rdev{};
        this->generator.seed(rdev());
    };

    RG_uniform_real(int seed, T min, T max){
    this->distribution = std::uniform_real_distribution<T>(min, max);
    this->generator.seed(seed);
    }

    std::vector<T> draw(size_t n){
        std::vector<T> data(n);
        std::generate(data.begin(), data.end(), [this]() { return distribution(generator); });
        return data;
    };

private:
    std::uniform_real_distribution<T> distribution;
    std::default_random_engine generator;
};

template <class T> class RG_uniform_int{
public:
    RG_uniform_int(){
        this->distribution = std::uniform_int_distribution<T>();
        std::random_device rdev{};
        this->generator.seed(rdev());
    };

    RG_uniform_int(int seed){
        this->distribution = std::uniform_int_distribution<T>();
        this->generator.seed(seed);
    };

    RG_uniform_int(T min, T max){
        this->distribution = std::uniform_int_distribution<T>(min, max);
        std::random_device rdev{};
        this->generator.seed(rdev());
    };

    RG_uniform_int(int seed, T min, T max){
        this->distribution = std::uniform_int_distribution<T>(min, max);
        this->generator.seed(seed);
    };

    std::vector<T> draw(size_t n) {
        std::vector<T> data(n);
        std::generate(data.begin(), data.end(), [this]() { return distribution(generator); });
        return data;
    }

private:
    std::uniform_int_distribution<T> distribution;
    std::default_random_engine generator;
};

#ifdef USE_CUDA
#define DEFAULT_SEED 7;
class RG_uniform_float{
public:
    RG_uniform_float(int devID){
//        printf("DEVICE: %i", devID);
//        printf("Allocating data for %i samples...\n", size);
        int seed = DEFAULT_SEED;
//        curandGenerator_t prngGPU;
        checkCudaErrors(curandCreateGenerator(&prngGPU, CURAND_RNG_PSEUDO_MTGP32));
        checkCudaErrors(curandSetPseudoRandomGeneratorSeed(prngGPU, seed));

//        printf("Seeding with %i ...\n", seed);
    }

    void alloc_mem(size_t size)
    {
        checkCudaErrors(cudaMalloc((void **)&d_Rand, size * sizeof(float)));
        h_RandGPU  = (float *)malloc(size * sizeof(float));
    }

    std::vector<float> draw(size_t size){

//        printf("Generating random numbers on GPU...\n\n");
        checkCudaErrors(curandGenerateUniform(prngGPU, (float *) d_Rand, size));

//        printf("\nReading back the results...\n");
        checkCudaErrors(cudaMemcpy(h_RandGPU, d_Rand, size * sizeof(float), cudaMemcpyDeviceToHost));

        std::vector<float> res;
        res.assign(h_RandGPU,h_RandGPU+size);
        return res;
    }

private:
    size_t size;
    curandGenerator_t prngGPU;
    float *d_Rand;
    float *h_RandGPU;

};

#endif

#endif //ECHELLESIMULATOR_RANDOM_GENERATOR_H
