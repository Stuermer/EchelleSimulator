//
// Created by stuermer on 3/27/17.
//

#ifndef ECHELLESIMULATOR_RANDOM_GENERATOR_H
#define ECHELLESIMULATOR_RANDOM_GENERATOR_H

#include <vector>
#include "stdlib.h"
#include <random>
#include <algorithm>
#include <type_traits>

template<class T>
class RandomGenerator {
public:
    virtual T draw() = 0;

    virtual std::vector<T> draw(size_t n) = 0;

    virtual T operator()() = 0;
};

template<class T>
class RG_uniform_real : RandomGenerator<T> {
public:
    RG_uniform_real() {
        this->distribution = std::uniform_real_distribution<T>();
        std::random_device rdev{};
        this->generator.seed(rdev());
    };


    RG_uniform_real(int seed) {
        this->distribution = std::uniform_real_distribution<T>();
        this->generator.seed(seed);
    };

    RG_uniform_real(T min, T max) {
        this->distribution = std::uniform_real_distribution<T>(min, max);
        std::random_device rdev{};
        this->generator.seed(rdev());
    };

    RG_uniform_real(int seed, T min, T max) {
        this->distribution = std::uniform_real_distribution<T>(min, max);
        this->generator.seed(seed);
    }

    std::vector<T> draw(size_t n) {
        std::vector<T> data(n);
        std::generate(data.begin(), data.end(), [this]() { return distribution(generator); });
        return data;
    };

private:
    std::uniform_real_distribution<T> distribution;
    std::default_random_engine generator;
};

template<class T>
class RG_uniform_int : RandomGenerator<T> {
public:
    RG_uniform_int() {
        this->distribution = std::uniform_int_distribution<T>();
        std::random_device rdev{};
        this->generator.seed(rdev());
    };

    RG_uniform_int(int seed) {
        this->distribution = std::uniform_int_distribution<T>();
        this->generator.seed(seed);
    };

    RG_uniform_int(T min, T max) {
        this->distribution = std::uniform_int_distribution<T>(min, max);
        std::random_device rdev{};
        this->generator.seed(rdev());
    };

    RG_uniform_int(int seed, T min, T max) {
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


template<class FloatType = double,
        class Generator = std::mt19937>
class piecewise_linear_RNG : public RandomGenerator<FloatType> {
public:
    typedef FloatType result_type;
    typedef Generator generator_type;
    typedef std::piecewise_constant_distribution<FloatType> distribution_type;

    // default constructor
    explicit piecewise_linear_RNG(std::vector<FloatType> values, std::vector<FloatType> weights,
                                  Generator &&_eng = Generator{std::random_device{}()}) : eng(std::move(_eng)),
                                                                                          dist(values.begin(),
                                                                                               values.end(),
                                                                                               weights.begin()) {}

    // construct from existing pre-defined engine
    explicit piecewise_linear_RNG(std::vector<FloatType> values, std::vector<FloatType> weights, const Generator &_eng)
            : eng(_eng), dist(values.begin(), values.end(), weights.begin()) {}

    // generate next random value in distribution (equivalent to next() in above code)
    result_type draw() { return dist(eng); }

    std::vector<result_type> draw(size_t n) {
        std::vector<result_type> data(n);
        std::generate(data.begin(), data.end(), [this]() { return dist(eng); });
        return data;
    }

    result_type operator()() { return dist(eng); }

private:
    generator_type eng;
    distribution_type dist;
};

template<class FloatType = double,
        class Generator = std::mt19937>
class discrete_RNG : public RandomGenerator<FloatType> {
public:
    typedef FloatType result_type;
    typedef Generator generator_type;
    typedef std::discrete_distribution<int> distribution_type;

    // default constructor
    explicit discrete_RNG(std::vector<FloatType> _values, std::vector<FloatType> _weights,
                          Generator &&_eng = Generator{std::random_device{}()}) : eng(std::move(_eng)), values(_values),
                                                                                  dist(_weights.begin(),
                                                                                       _weights.end()) {}

    // construct from existing pre-defined engine
    explicit discrete_RNG(std::vector<FloatType> _values, std::vector<FloatType> _weights, const Generator &_eng)
            : eng(_eng), values(_values), dist(_weights.begin(), _weights.end()) {}

    // generate next random value in distribution (equivalent to next() in above code)
    result_type draw() { return values[dist(eng)]; }

    std::vector<result_type> draw(size_t n) {
        std::vector<result_type> data(n);
        std::generate(data.begin(), data.end(), [this]() { return values[dist(eng)]; });
        return data;
    }

    result_type operator()() { return values[dist(eng)]; }

private:
    generator_type eng;
    std::vector<result_type> values;
    distribution_type dist;
};

#endif //ECHELLESIMULATOR_RANDOM_GENERATOR_H
