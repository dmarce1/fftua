#include "fftu.hpp"
#include "util.hpp"
#include <unordered_map>
#include <cstring>


void fft_raders_indices_real(int* I, int N) {
	std::vector<int> J(N);
	const auto& gq = raders_gq(N);
	for (int n = 0; n < N - 1; n++) {
		J[n] = I[gq[n]];
	}
	J[N - 1] = I[0];
	std::memcpy(I, J.data(), N * sizeof(int));
	fft_indices_real(I, N - 1);
}

void fft_raders_dht(fft_simd4* X, int N) {
	const int M = N - 1;
	workspace<fft_simd4> ws;
	auto Y = ws.create(M);
	const auto& gq = raders_gq(N);
	const auto& ginvq = raders_ginvq(N);
	const auto& W = raders_twiddle_real(N, M);
	fft_simd4 xo = X[N - 1];
	fft_simd4 x0 = xo;
	fht(X, M);
	x0 += X[0];
	Y[0] = X[0] * W[0];
	if (M % 2 == 0) {
		Y[M / 2] = X[M / 2] * W[M / 2];
	}
	for (int q = 1; q < M - q; q++) {
		const auto yp = (X[q] + X[M - q]) * 0.5;
		const auto ym = (X[q] - X[M - q]) * 0.5;
		const auto pos = yp * W[q] + ym * W[M - q];
		const auto neg = yp * W[M - q] - ym * W[q];
		Y[q] = pos;
		Y[M - q] = neg;
	}
	fft_scramble_real(Y.data(), M);
	fht(Y.data(), M);
	X[0] = x0;
	for (int p = 0; p < N - 1; p++) {
		X[ginvq[p]] = Y[p] + xo;
	}
	ws.destroy(std::move(Y));
}

void fft_raders_dht(double* X, int N, bool padded) {
	const int M = padded ? compute_padding(N - 1) : N - 1;
	const int ysize = round_up(M, SIMD_SIZE);
	workspace<double> ws;
	auto Y = ws.create(ysize);
	const auto& gq = raders_gq(N);
	const auto& ginvq = raders_ginvq(N);
	const auto& W = raders_twiddle_real(N, M);
	double xo = X[0];
	double x0 = xo;
	for (int n = 0; n < N - 1; n++) {
		Y[n] = X[gq[n]];
	}
	for (int n = N - 1; n < M; n++) {
		Y[n] = 0.0;
	}
	fht(Y.data(), M);
	x0 += Y[0];
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
	ws.destroy(std::move(Y));
}

void fft_raders_real(double* X, int N, bool padded) {
	fft_raders_dht(X, N, padded);
	fht2fft(X, N);
}

void fft_raders_real(fft_simd4* X, int N) {
	fft_raders_dht(X, N);
	fht2fft(X, N);
}
