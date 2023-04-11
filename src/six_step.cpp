#include "fftu.hpp"
#include "util.hpp"
#include <unordered_map>
#include <cmath>
#include <vector>
#include <stack>
#include <cstring>

void fft_six_step(complex<fft_simd4>* X, int N) {
	const int M = (1 << (ilogb(N) >> 1));
	const auto& W = twiddles(N);
	for (int n = 0; n < M; n++) {
		for (int m = n + 1; m < M; m++) {
			std::swap(X[M * n + m], X[M * m + n]);
		}
	}
	for (int n = 0; n < M; n++) {
		fft_scramble(X + n * M, M);
		fft(X + n * M, M);
	}
	for (int n = 1; n < M; n++) {
		const auto w0 = W[n];
		auto w = w0;
		for (int m = 1; m < M; m++) {
			X[M * n + m] *= w;
			w *= w0;
		}
	}
	for (int n = 0; n < M; n++) {
		for (int m = n + 1; m < M; m++) {
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
