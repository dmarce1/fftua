#include "fftu.hpp"
#include "util.hpp"
#include <unordered_map>

template<class T>
void fht(T* X, int N) {
	fft_real(X, N);
	fft2dht(X, N);
}

void fft_raders_dht(double* X, int N, bool padded) {
	static thread_local std::unordered_map<int, std::vector<double>> cache;
	std::vector<double>& Y = cache[N];
	const int M = padded ? compute_padding(N - 1) : N - 1;
	const int ysize = round_up(M, SIMD_SIZE);
	Y.resize(ysize);
	const auto& gq = raders_gq(N);
	const auto& ginvq = raders_ginvq(N);
	const auto& W = raders_twiddle_real(N, M);
	double xo = X[0];
	double x0 = xo;
	for (int n = 0; n < N - 1; n++) {
		Y[n] = X[gq[n]];
		x0 += X[gq[n]];
	}
	for (int n = N - 1; n < M; n++) {
		Y[n] = 0.0;
	}
	fht(Y.data(), M);
	Y[0] *= W[0];
	if (M % 2 == 0) {
		Y[M / 2] *= W[M / 2];
	}
	for (int q = 1; q < M - q; q++) {
		const auto yp = (Y[q] + Y[M - q]) * 0.5;
		const auto ym = (Y[q] - Y[M - q]) * 0.5;
		const auto pos = W[q] * yp + W[M - q] * ym;
		const auto neg = W[M - q] * yp - W[q] * ym;
		Y[q] = pos;
		Y[M - q] = neg;
	}
	fht(Y.data(), M);
	X[0] = x0;
	for (int p = 0; p < N - 1; p++) {
		X[ginvq[p]] = Y[p] + xo;
	}
}

void fft_raders_real(double* X, int N, bool padded) {
	fft_raders_dht(X, N, padded);
	dht2fft(X, N);
}
