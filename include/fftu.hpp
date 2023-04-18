/*
 * fftu.hpp
 *
 *  Created on: Apr 10, 2023
 *      Author: dmarce1
 */

#ifndef FFTU_HPP_
#define FFTU_HPP_

#include "types.hpp"
#include <vector>

#include "sfft.hpp"

void fft_cooley_tukey(int N1, complex<fft_simd4>* X, int N);
void fft_cooley_tukey(int N1, complex<double>* X, int N);
void fft_cooley_tukey_indices(int N1, int* I, int N);
void fft_conjugate(int N1, complex<fft_simd4>* X, int N);
void fft_conjugate(int N1, complex<double>* X, int N);
void fft_conjugate_indices(int N1, int* I, int N);
void fft_split(int R, complex<double>* X, int N);
void fft_split(int R, complex<fft_simd4>* X, int N);
void fft_split_indices(int R, int* I, int N);
void fft_split_conjugate(int R, complex<double>* X, int N);
void fft_split_conjugate(int R, complex<fft_simd4>* X, int N);
void fft_split_conjugate_indices(int R, int* I, int N);
const std::vector<int>& fft_indices(int N);
const std::vector<int>& fft_inv_indices(int N);
void fft(complex<fft_simd4>* X, int N);
void fft(complex<double>* X, int N);
void fft_indices(int*, int);
void fft_scramble(complex<fft_simd4>* X, int N);
void fft_scramble(complex<double>* X, int N);
void fft_six_step_indices(int* I, int N);
void fft_six_step(complex<double>* X, int N);
void fft_six_step(complex<fft_simd4>* X, int N);
void fft_raders_indices(int* I, int N);
void fft_raders(complex<fft_simd4>* X, int N, bool scramble = false);
void fft_raders(complex<double>* X, int N, bool scramble = false);
void fft_permute(const std::vector<int>&, complex<double>* X);
void fft_permute(const std::vector<int>&, complex<fft_simd4>* X);

inline double rand1() {
	return (rand() + 0.5) / (RAND_MAX + 1.0);
}

#endif /* FFTU_HPP_ */
