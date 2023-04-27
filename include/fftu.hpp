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
#include <string>

#include "sfft.hpp"

void fft_twoforone_real(fft_simd4* X, int N);
void fft_twoforone_real(double* X, int N);
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
void fft_scramble_inv(complex<fft_simd4>* X, int N);
void fft_scramble_inv(complex<double>* X, int N);
void fft_six_step_indices(int* I, int N);
void fft_six_step(complex<double>* X, int N);
void fft_six_step(complex<fft_simd4>* X, int N);
void fft_raders_indices(int* I, int N);
void fft_raders_padded_indices(int* I, int N);
void fft_raders(complex<fft_simd4>* X, int N);
void fft_raders(complex<double>* X, int N, bool padded);
void fft_permute(const std::vector<int>&, complex<double>* X);
void fft_permute(const std::vector<int>&, complex<fft_simd4>* X);
void fft_raders_padded(complex<fft_simd4>* X, int N);
void fft_bluestein(complex<double>* X, int N);
void fft_real(double* X, int N);
void fft_split_real_indices(int R, int* I, int N);
void fft_split_real(int R, double* X, int N);
void fft_split_real(int R, fft_simd4* X, int N);
void fft_twoforone_real1(fft_simd4* X, int N);
void fft_fourfortwo_real(double* X, int N);

#define FFT_SPLIT 0
#define FFT_SPLIT_CONJ 4
#define FFT_CT 1
#define FFT_CONJ 2
#define FFT_6 3
#define FFT_RADERS 5
#define FFT_RADERS_PADDED 6
#define FFT_241 7
#define FFT_PFAC 8
#define FFT_BLUE 9
#define FFT_442 10

struct fft_method {
	int type;
	int R;
};

std::string fft_method_string(const fft_method method);

inline double rand1() {
	return (rand() + 0.5) / (RAND_MAX + 1.0);
}

void fht_radix2(double* x, int N, bool scramble = true);

template<class T>
void fft_permute1(const std::vector<int>& I, T* X) {
	const int N = I.size();
	static std::vector<bool> flag;
	flag.resize(N, false);
	for (int n = 0; n < N; n++) {
		if (!flag[n]) {
			flag[n] = true;
			T tmp = X[n];
			int m = I[n];
			while (m != n) {
				flag[m] = true;
				std::swap(tmp, X[m]);
				m = I[m];
			}
			X[n] = tmp;
		}
	}
	flag.resize(0);

}


inline void fft_permute(const std::vector<int>& I, complex<double>* X) {
	fft_permute1(I, X);
}

inline void fft_permute(const std::vector<int>& I, complex<fft_simd4>* X) {
	fft_permute1(I, X);
}

inline void fft_permute(const std::vector<int>& I, double* X) {
	fft_permute1(I, X);
}

inline void fft_permute(const std::vector<int>& I, fft_simd4* X) {
	fft_permute1(I, X);
}


void fft_cooley_tukey_real(int N1, complex<double>* X, int N);
void fft_cooley_tukey_real(int N1, double* X, int N);
void fft_cooley_tukey_indices_real(int N1, int* I, int N);
void fft_conjugate_real(int N1, fft_simd4* X, int N);
void fft_split_conjugate_real(int N1, fft_simd4* X, int N);
void fft_conjugate_indices_real(int N1, int* I, int N);
void fft_split_conjugate_indices_real(int N1, int* I, int N);
void fft_cooley_tukey_real(int N1, fft_simd4* X, int N);
void fft_cooley_tukey_real(int N1, double* X, int N);
void fft_indices_real(int* I, int N);
void fft_real(fft_simd4* X, int N);
void fft_real(double* X, int N);
void fft_raders_real(double* X, int N, bool padded);
const std::vector<int>& fft_indices_real(int N);
void fft_raders_real(fft_simd4* X, int N);
void fft_scramble_real(fft_simd4* X, int N);
void fft_scramble_real(double* X, int N);
void fft_raders_prime_factor_real(int N1, double* X, int N);

template<class T>
void fft2fht(T* x, int N) {
	for (int n = 1; n < N - n; n++) {
		const auto r = x[n];
		const auto i = x[N - n];
		x[n] = r - i;
		x[N - n] = r + i;
	}
}

template<class T>
void fht2fft(T* x, int N) {
	for (int j = 1; j < N - j; j++) {
		const auto p = x[j];
		const auto n = x[N - j];
		x[j] = T(0.5) * (p + n);
		x[N - j] = T(0.5) * (n - p);
	}
}


template<class T>
void fht(T* X, int N) {
	fft_real(X, N);
	fft2fht(X, N);
}

#endif /* FFTU_HPP_ */
