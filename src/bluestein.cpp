#include "fftu.hpp"
#include "util.hpp"
#include <stack>
#include <cstring>

void fft_bluestein(complex<double>* X, int N) {
	static std::stack<std::vector<complex<double>>>astack;
	std::vector<complex<double>> a;
	if( astack.size()) {
		a = std::move(astack.top());
		astack.pop();
	}
	const int M = compute_padding(N);
	const auto m = bluestein_multiplier(N);
	const auto h = bluestein_filter(N, M);
	a.resize(M);
	const int end = (N % 2 == 1) ? N - 1 : N;
	for (int n = 0; n < end; n += 2) {
		a[n] = X[n];
		a[n + 1] = X[n + 1];
		*((__m256d *) &a[n]) = mul(*((__m256d *) &a[n]), *((__m256d *) &m[n]));
	}
	if (N % 2 == 1) {
		a[N - 1] = X[N - 1] * m[N - 1];
	}
	for (int n = N; n < M; n++) {
		a[n] = 0.0;
	}
	fft(a.data(), M);
	for (int k = 0; k < M; k += 2) {
		*((__m256d *) &a[k]) = mul(*((__m256d *) &a[k]), *((__m256d *) &h[k]));
	}
	for (int n = 1; n < M - n; n++) {
		std::swap(a[n], a[M - n]);
	}
	fft(a.data(), M);
	for (int n = 0; n < end; n += 2) {
		*((__m256d *) &a[n]) = mul(*((__m256d *) &a[n]), *((__m256d *) &m[n]));
		X[n] = a[n];
		X[n + 1] = a[n + 1];
	}
	if (N % 2 == 1) {
		X[N - 1] = a[N - 1] * m[N - 1];
	}
	astack.push(std::move(a));
}
