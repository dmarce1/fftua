#include "fftu.hpp"
#include "util.hpp"
#include <cstring>

void fft_raders_indices(int* I, int N) {
	std::vector<int> J(N);
	const auto& gq = raders_gq(N);
	for (int n = 0; n < N - 1; n++) {
		J[n] = I[gq[n]];
	}
	J[N - 1] = I[0];
	std::memcpy(I, J.data(), N * sizeof(int));
	fft_indices(I, N - 1);
}

void fft_raders(complex<fft_simd4>* X, int N, bool scramble) {
	static thread_local std::unordered_map<int, std::vector<complex<fft_simd4>>>cache;
	std::vector<complex<fft_simd4>>& Y = cache[N];
	Y.resize(N - 1);
	const auto& gq = raders_gq(N);
	const auto& tw = raders_twiddle(N);
	complex<fft_simd4> xo = scramble ? X[0] : X[N - 1];
	complex<fft_simd4> x0 = xo;
	if( scramble ) {
		for (int n = 0; n < N - 1; n++) {
			Y[n] = X[gq[n]];
		}
		fft_scramble(Y.data(), N - 1);
		fft(Y.data(), N - 1);
		x0 += Y[0];
		for (int q = 0; q < N - 1; q++) {
			Y[q] *= tw[q];
		}
		fft_scramble(Y.data(), N - 1);
	} else {
		fft(X, N - 1);
		x0 += X[0];
		for (int q = 0; q < N - 1; q++) {
			X[q] *= tw[q];
		}
		const auto& I = fft_inv_indices(N - 1);
		for (int q = 0; q < N - 1; q++) {
			Y[I[q]] = X[q];
		}
	}
	fft(Y.data(), N - 1);
	for (int p = 0; p < N - 1; p++) {
		Y[p] += xo;
	}
	X[0] = x0;
	for (int p = 0; p < N - 1; p++) {
		X[gq[p]] = Y[p];
	}
}

void fft_raders(complex<double>* X, int N, bool scramble) {
	static thread_local std::unordered_map<int, std::vector<complex<double>>>cache;
	std::vector<complex<double>>& Y = cache[N];
	const int ysize = round_up(N - 1, SIMD_SIZE);
	Y.resize(ysize);
	fft_simd4* Z = (fft_simd4*) Y.data();
	const auto& gq = raders_gq(N);
	const auto& tw = raders_twiddle(N);
	complex<double> xo = X[0];
	complex<double> x0 = xo;
	fft_simd4 z0;
	for( int n = 0; n < SIMD_SIZE; n += 2) {
		z0[n] = xo.real();
	}
	for( int n = 1; n < SIMD_SIZE; n += 2) {
		z0[n] = xo.imag();
	}
	for (int n = 0; n < N - 1; n++) {
		Y[n] = X[gq[n]];
	}
	fft(Y.data(), N - 1);
	x0 += Y[0];
	for (int q = 0; q < (N - 1); q++) {
		Y[q] *= tw[q];
	}
	fft(Y.data(), N - 1);
	for (int p = 0; p < 2 * ysize / SIMD_SIZE; p++) {
		Z[p] += z0;
	}
	X[0] = x0;
	for (int p = 0; p < N - 1; p++) {
		X[gq[p]] = Y[p];
	}
}

