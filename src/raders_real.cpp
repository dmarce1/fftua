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

double fht_convolve(double* x, const double* h, int N) {
	fht(x, N);
	double rc = x[0];
	x[0] *= h[0];
	if (N % 2 == 0) {
		x[N / 2] *= h[N / 2];
	}
	for (int q = 1; q < N - q; q++) {
		const auto yp = (x[q] + x[N - q]) * 0.5;
		const auto ym = (x[q] - x[N - q]) * 0.5;
		x[q] = h[q] * yp + h[N - q] * ym;
		x[N - q] = h[N - q] * yp - h[q] * ym;
	}
	fht(x, N);
	return rc;
}

void fft_raders_dht(double* X, int N, bool padded) {
	auto prime_fac = prime_factorization(N);
	assert(prime_fac.size() == 1);
	int P = prime_fac.begin()->first;
	int c = prime_fac.begin()->second;
	int L = std::pow(P, c - 1);
	int M = L * (P - 1);
	int K = padded ? compute_padding(M) : M;
	const auto& W = raders_twiddle_real(N, K);
	const auto& gq = raders_gq(N);
	const auto& ginvq = raders_ginvq(N);
	workspace<double> ws;
	auto x2 = ws.create(K);
	if (L > 1) {
		auto x0 = ws.create(L);
		auto x1 = ws.create(L);
		for (int n2 = 0; n2 < L; n2++) {
			x0[n2] = 0.0;
			for (int n1 = 0; n1 < P; n1++) {
				x0[n2] += X[L * n1 + n2];
			}
		}
		for (int n1 = 0; n1 < L; n1++) {
			x1[n1] = X[P * n1];
		}
		for (int q = 0; q < M; q++) {
			x2[q] = X[gq[q]];
		}
		for (int q = M; q < K; q++) {
			x2[q] = 0.0;
		}
		fht(x0.data(), L);
		fht(x1.data(), L);
		fht_convolve(x2.data(), W.data(), M);
		for (int k = 0; k < N; k++) {
			X[k] = 0.0;
		}
		for (int k1 = 0; k1 < L; k1++) {
			X[P * k1] += x0[k1];
		}
		for (int k1 = 0; k1 < P; k1++) {
			for (int k2 = 0; k2 < L; k2++) {
				if (k2 % P != 0) {
					X[L * k1 + k2] += x1[k2];
				}
			}
		}
		ws.destroy(std::move(x0));
		ws.destroy(std::move(x1));
		for (int p = 0; p < M; p++) {
			X[ginvq[p]] += x2[p];
		}
	} else {
		for (int q = 0; q < M; q++) {
			x2[q] = X[gq[q]];
		}
		for (int q = M; q < K; q++) {
			x2[q] = 0.0;
		}
		auto xo = X[0];
		X[0] = xo + fht_convolve(x2.data(), W.data(), K);
		for (int p = 0; p < M; p++) {
			X[ginvq[p]] = xo + x2[p];
		}
	}
	ws.destroy(std::move(x2));
}

void fft_raders_real(double* X, int N, bool padded) {
	fft_raders_dht(X, N, padded);
	fht2fft(X, N);
}

void fft_raders_real(fft_simd4* X, int N) {
	fft_raders_dht(X, N);
	fht2fft(X, N);
}
