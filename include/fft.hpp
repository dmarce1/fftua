#pragma once

#include "sfft.hpp"
#include "types.hpp"
#include <stdio.h>
#include <vector>

#define FFT_COOLEY_TUKEY 0

struct fft_method_t {
	int type;
	int radix;
};

void fft_cooley_tukey(int N1, complex<double>* X, int N);
void fft_cooley_tukey(int N1, complex<fft_simd4>* X, int N);

const std::vector<int>& fft_indices(int N);
void fft_cooley_tukey_indices(int N1, int* I, int N);

fft_method_t fft_select(int N);

template<class T>
void fft_scramble(complex<T>* X, int N) {
	const auto I = fft_indices(N);
	static thread_local std::vector<bool> flag;
	flag.resize(N, false);
	for (int n = 0; n < N; n++) {
		if (!flag[n]) {
			flag[n] = true;
			auto tmp = X[n];
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

template<class T>
void fft(complex<T>* X, int N, bool scram = true) {
	const auto method = fft_select(N);
	switch (method.type) {
	case FFT_COOLEY_TUKEY:
		if (scram) {
			fft_scramble(X, N);
		}
		fft_cooley_tukey(method.radix, X, N);
		break;
	default:
		assert(false);
		abort();
	}
}

void fft_indices(int* I, int N);
