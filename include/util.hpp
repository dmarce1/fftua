/*
 * util.hpp
 *
 *  Created on: Mar 5, 2023
 *      Author: dmarce1
 */

#ifndef UTIL_HPP_
#define UTIL_HPP_

#include <vector>
#include <set>
#include <unordered_map>
#include <map>
#include <functional>
#include <vector>
#include <complex>
#include "types.hpp"
#include <chrono>

#define SIMD_SIZE 4

class timer {
	std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
	double time;
public:
	inline timer() {
		time = 0.0;
	}
	inline void stop() {
		std::chrono::time_point<std::chrono::high_resolution_clock> stop_time = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> dur = stop_time - start_time;
		time += dur.count();
	}
	inline void start() {
		start_time = std::chrono::high_resolution_clock::now();
	}
	inline void reset() {
		time = 0.0;
	}
	inline double read() {
		return time;
	}
};

inline bool close2(double a, double b) {
	return std::abs(a - b) < 1.0e-10;
}
int least_prime_factor(int N);
int mod(int a, int b);
bool power_of(int N, int M);
std::vector<std::vector<int>> nchoosek(int n, int k);
bool is_prime(int n);
int greatest_prime_factor(int N);
std::map<int, int> prime_factorization(int N);
int mod_pow(int a, int b, int m);
int mod_inv(int a, int m);
int generator(int N);
std::vector<int> raders_ginvq(int N);
const std::vector<int> raders_gq(int N);
double fftw(std::vector<std::complex<double>>& x);
const std::vector<complex<double>> raders_twiddle(int N, int M, bool padded);
bool are_coprime(int a, int b);
std::vector<complex<double>> chirp_z_filter(int N);
const std::vector<complex<double>>& twiddles(int N);
inline int round_down(int i, int m) {
	return m * (i / m);
}
inline int round_up(int i, int m) {
	return m * (((i - 1) / m) + 1);
}
const std::vector<std::vector<complex<fft_simd4>>>& vector_twiddles(int N1, int N2);
void destroy_scratch(std::vector<complex<fft_simd4>>&& space);
std::vector<complex<fft_simd4>> create_scratch(int N);

#endif
