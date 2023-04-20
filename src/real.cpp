#include "fftu.hpp"
#include "sfft.hpp"
#include "util.hpp"
#include <cstring>


void fft_real(int N1, double* X, int N) {
	if (N <= SFFT_NMAX) {
		sfft_real(X, N);
		return;
	}
	const int N2 = N / N1;
	const auto& W = twiddles(N);
	std::vector<double> Y;
	complex<double> p[N1];
	Y.resize(N);

	for (int n1 = 0; n1 < N1; n1++) {
		for (int n2 = 0; n2 < N2; n2++) {
			Y[N2 * n1 + n2] = X[N1 * n2 + n1];
		}
	}

	for (int n1 = 0; n1 < N1; n1++) {
		fft_real(N1, &Y[n1 * N2], N2);
	}

	for (int k2 = 1; k2 < (N2 + 1) / 2; k2++) {
		for (int n1 = 1; n1 < N1; n1++) {
			const auto& w = W[n1 * k2];
			auto& x = Y[N2 * n1 + mod(+k2, N2)];
			auto& y = Y[N2 * n1 + mod(-k2, N2)];
			auto tmp = x;
			x = x * w.real() - y * w.imag();
			y = tmp * w.imag() + y * w.real();
		}
	}

	std::memcpy(X, Y.data(), sizeof(double) * N);
	{
		double q[N1];
		for (int n1 = 0; n1 < N1; n1++) {
			q[n1] = Y[N2 * n1];
		}
		sfft_real(q, N1);
		X[0] = q[0];
		for (int k1 = 1; k1 < (N1 + 1) / 2; k1++) {
			X[0 + N2 * k1] = q[mod(+k1, N1)];
			X[N - N2 * k1] = q[mod(-k1, N1)];
		}
		if (N1 % 2 == 0) {
			X[N / 2] = q[N1 / 2];
		}
	}
	for (int k2 = 1; k2 < (N2 + 1) / 2; k2++) {
		for (int n1 = 0; n1 < N1; n1++) {
			p[n1].real() = Y[N2 * n1 + mod(+k2, N2)];
			p[n1].imag() = Y[N2 * n1 + mod(-k2, N2)];
		}
		sfft_complex((double*) p, N1);
		for (int k1 = 0; k1 < (N1 + 1) / 2; k1++) {
			const int k = N2 * k1 + k2;
			X[k] = p[k1].real();
			X[N - k] = p[k1].imag();
		}
		for (int k1 = (N1 + 1) / 2; k1 < N1; k1++) {
			const int k = N2 * k1 + k2;
			X[N - k] = p[k1].real();
			X[k] = -p[k1].imag();
		}
	}
	if (N2 % 2 == 0) {
		double q[N1];
		for (int n1 = 0; n1 < N1; n1++) {
			q[n1] = Y[N2 * n1 + N2 / 2];
		}
		sfft_skew(q, N1);
		for (int k1 = 0; k1 < (N1 + 1) / 2; k1++) {
			X[0 + (N2 * k1 + N2 / 2)] = q[k1];
			X[N - (N2 * k1 + N2 / 2)] = q[N1 - k1 - 1];
		}
		if (N1 % 2 == 1) {
			X[N / 2] = q[N1 / 2];
		}
	}
}
