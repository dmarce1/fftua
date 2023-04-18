#include "fftu.hpp"
#include "util.hpp"
#include <cstring>

template<class T>
void fft_raders1(complex<T>* X, int N, bool scramble) {
	static thread_local std::unordered_map<int, std::vector<complex<T>>>cache;
	std::vector<complex<T>>& Y = cache[N];
	Y.resize(N - 1);
	const auto& gq = raders_gq(N);
	const auto& ginvq = raders_ginvq(N);
	const auto& tw = raders_twiddle(N);
	constexpr bool simd = !std::is_same<double, T>::value;
	complex<T> xo = (simd && !scramble) ? X[N - 1] : X[0];
	complex<T> x0 = xo;
	if( simd ) {
		if( scramble ) {
			for (int n = 1; n < N; n++) {
				x0 += X[n];
			}
			for (int n = 0; n < N - 1; n++) {
				Y[n] = X[gq[n]];
			}
			fft_scramble(Y.data(), N - 1);
			fft(Y.data(), N - 1);
			for (int q = 0; q < N - 1; q++) {
				Y[q] *= tw[q];
			}
			for (int q = 1; q < N - 1 - q; q++) {
				std::swap(Y[q], Y[N - 1 - q]);
			}
			fft_scramble(Y.data(), N - 1);
		} else {
			for (int n = 0; n < N - 1; n++) {
				x0 += X[n];
			}
			fft(X, N - 1);
			for (int q = 0; q < N - 1; q++) {
				X[q] *= tw[q];
			}
			const auto& I = fft_inv_indices(N - 1);
			Y[I[0]] = X[0];
			for (int q = 1; q < N - 1; q++) {
				Y[I[q]] = X[N - 1 - q];
			}
		}
	} else {
		for (int n = 1; n < N; n++) {
			x0 += X[n];
		}
		for (int n = 0; n < N - 1; n++) {
			Y[n] = X[gq[n]];
		}
		fft(Y.data(), N - 1);
		for (int q = 0; q < N - 1; q++) {
			Y[q] *= tw[q];
		}
		for (int q = 1; q < N - 1 - q; q++) {
			std::swap(Y[q], Y[N - 1 - q]);
		}
	}
	fft(Y.data(), N - 1);
	for (int p = 0; p < N - 1; p++) {
		Y[p] += xo;
	}
	X[0] = x0;
	for (int p = 0; p < N - 1; p++) {
		X[ginvq[p]] = Y[p];
	}
}

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
	fft_raders1<fft_simd4>(X, N, scramble);
}

void fft_raders(complex<double>* X, int N, bool scramble) {
	fft_raders1<double>(X, N, false);
}

