#pragma once

#include "types.hpp"
#include <map>
#include <vector>

bool is_prime(int);
const std::map<int, int>& prime_factorization(int);
const std::vector<complex<double>>& twiddles(int);
double fftw(std::vector<complex<double>>&);
double rand1();
int greatest_prime_factor(int);
