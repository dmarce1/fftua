#include "fftu.hpp"
#include "util.hpp"
#include <unordered_map>
#include <cmath>
#include <vector>
#include <stack>
#include <cstring>

void fft_six_step_indices(int* I, int N) {
	const int M = lround(sqrt(N));
	std::vector<int> J(N);
	for (int n = 0; n < M; n++) {
		for (int m = 0; m < M; m++) {
			J[M * n + m] = I[M * m + n];
		}
	}
	const auto K = fft_indices(M);
	for (int n = 0; n < M; n++) {
		for (int m = 0; m < M; m++) {
			I[M * n + m] = J[M * n + K[m]];
		}
	}
}

template<class T>
void fft_six_step1(complex<T>* X, int N) {
	const int M = lround(sqrtl(N));
	const auto& W = twiddles(N);
	for (int n = 0; n < M; n++) {
		fft(X + n * M, M);
	}
	for (int n = 0; n < M; n++) {
		X[M * n + n] *= W[n * n];
		for (int m = n + 1; m < M; m++) {
			const auto Mnm = M * n + m;
			const auto Mmn = M * m + n;
			const auto& w = W[n * m];
			X[Mnm] *= w;
			X[Mmn] *= w;
			std::swap(X[Mnm], X[Mmn]);
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

void fft_six_step(complex<double>* X, int N) {
	fft_six_step1(X, N);
}

void fft_six_step(complex<fft_simd4>* X, int N) {
	fft_six_step1(X, N);
}
