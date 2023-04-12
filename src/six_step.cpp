#include "fftu.hpp"
#include "util.hpp"
#include <unordered_map>
#include <cmath>
#include <vector>
#include <stack>
#include <cstring>

void fft_six_step_indices(int* I, int N) {
	const int M = (1 << (ilogb(N) >> 1));
	std::vector<int> J(N);
	for (int n = 0; n < M; n++) {
		for (int m = 0; m < M; m++) {
			J[M * n + m] = I[M * m + n];
		}
	}
	const auto K = fft_indices(M);
	for (int n = 0; n < M; n++) {
		for (int m = 0; m < M; m++) {
			I[M * n + K[m]] = J[M * n + m];
		}
	}
}

void fft_six_step(complex<fft_simd4>* X, int N) {
	const int M = (1 << (ilogb(N) >> 1));
	const auto& W = twiddles(N);
	for (int n = 0; n < M; n++) {
		fft(X + n * M, M);
	}
	for (int n = 1; n < M; n++) {
		X[M * n + n] *= W[n * n];
	}
	for (int n = 0; n < M; n++) {
		for (int m = n + 1; m < M; m++) {
			X[M * n + m] *= W[n * m];
			X[M * m + n] *= W[n * m];
			std::swap(X[M * n + m], X[M * m + n]);
		}
	}
	for (int n = 0; n < M; n++) {
		fft_scramble(X + n * M, M);
		fft(X + n * M, M);
	}
	for (int n = 0; n < M; n++) {
		for (int m = n + 1; m < M; m++) {
			std::swap(X[M * n + m], X[M * m + n]);
		}
	}
}
