
#include "catch.hpp"
#include "../serial_src/rng.h"
#include <vector>
#include <numeric>
#include <cmath>
#include<iostream>

TEST_CASE("RNG reproducibility with same seed", "[RNG]") {
    RNG rng1(42);
    RNG rng2(42);

    std::vector<double> seq1(10);
    std::vector<double> seq2(10);

    for (int i = 0; i < 10; ++i) {
        seq1[i] = rng1.uniform01();
        seq2[i] = rng2.uniform01();
    }

    REQUIRE(seq1 == seq2);
}

TEST_CASE("RNG different seeds produce different first values", "[RNG]") {
    RNG rng1(1);
    RNG rng2(2);

    double v1 = rng1.uniform01();
    double v2 = rng2.uniform01();

    REQUIRE(v1 != Approx(v2));
}

TEST_CASE("RNG uniform01 range", "[RNG]") {
    RNG rng(123);

    for (int i = 0; i < 1000; ++i) {
        double v = rng.uniform01();
        REQUIRE(v >= 0.0);
        REQUIRE(v < 1.0);
    }
}

TEST_CASE("RNG uniformInt range", "[RNG]") {
    RNG rng(456);
    int low = 5;
    int high = 10;

    for (int i = 0; i < 1000; ++i) {
        int v = rng.uniformInt(low, high);
        REQUIRE(v >= low);
        REQUIRE(v <= high);
    }
}

TEST_CASE("RNG uniform(a,b) reproducibility", "[RNG]") {
    RNG rng1(7);
    RNG rng2(7);
    double a = 2.0;
    double b = 5.0;

    std::vector<double> seq1(10);
    std::vector<double> seq2(10);

    for (int i = 0; i < 10; ++i) {
        seq1[i] = rng1.uniform(a, b);
        seq2[i] = rng2.uniform(a, b);
    }

    REQUIRE(seq1 == seq2);
}

TEST_CASE("RNG adjacent correlation is low", "[RNG][correlation]") {
    RNG rng(909);
    const int N = 1000000;
    std::vector<double> seq(N);
    for (int i = 0; i < N; ++i) {
        seq[i] = rng.uniform01();
    }
    // Compute means
    double sumX = std::accumulate(seq.begin(), seq.end() - 1, 0.0);
    double sumY = std::accumulate(seq.begin() + 1, seq.end(), 0.0);
    double meanX = sumX / (N - 1);
    double meanY = sumY / (N - 1);
    
    // Compute variances and covariance
    double varX = 0.0, varY = 0.0, cov = 0.0;
    for (int i = 0; i < N - 1; ++i) {
        double dx = seq[i] - meanX;
        double dy = seq[i + 1] - meanY;
        varX += dx * dx;
        varY += dy * dy;
        cov += dx * dy;
    }
    varX /= (N - 1);
    varY /= (N - 1);
    cov /= (N - 1);
    double corr = cov / std::sqrt(varX * varY);

    // For a good RNG, correlation at lag 1 should be near zero
    REQUIRE(std::fabs(corr) < 0.005);
}


TEST_CASE("RNG uniform01 statistical uniformity", "[RNG][uniformity]") {
    RNG rng(205);
    const int N = 1000000;
    const int bins = 10;
    std::vector<int> counts(bins, 0);
    // Generate samples and bin them
    for (int i = 0; i < N; ++i) {
        double v = rng.uniform01();
        int bin = static_cast<int>(v * bins);
        if (bin == bins) bin = bins - 1; // edge case v==1.0
        counts[bin]++;
    }
    // Expected count per bin
    double expected = static_cast<double>(N) / bins;
    // Compute chi-squared statistic
    double chi2 = 0.0;
    for (int i = 0; i < bins; ++i) {
        double diff = counts[i] - expected;
        chi2 += diff * diff / expected;
    }
    // Output chi-squared value for reference
    std::cout<<"statistical uniformity test using a chi-squared goodness-of-fit -> Chi² = " << chi2<<std::endl;
    // Degrees of freedom = bins - 1 = 9; critical value at p=0.05 is ~16.919
    double critical = 16.919;
    REQUIRE(chi2 < critical);
}
