#include "fftu.hpp"
#include "util.hpp"

void convolve(double* x, const double* h, int N) {
	fht(x, N);
	x[0] *= h[0];
	if (N % 2 == 0) {
		x[N / 2] *= h[N / 2];
	}
	const double c0 = 0.5 / N;
	for (int q = 1; q < N - q; q++) {
		const auto yp = (x[q] + x[N - q]) * c0;
		const auto ym = (x[q] - x[N - q]) * c0;
		const auto pos = yp * h[q] + ym * h[N - q];
		const auto neg = yp * h[N - q] - ym * h[q];
		x[q] = pos;
		x[N - q] = neg;
	}
	fht(x, N);
}

void fft_raders_prime_factor_real(int N1, double* X, int N) {
	int N2 = N / N1;
	if (N2 % 2 == 0) {
		std::swap(N1, N2);
	}
	workspace<double> ws;
	workspace<complex<double>> cws;
	auto Y = ws.create(N);
	auto z = cws.create(N1);
	for (int n = 0; n < N; n++) {
		const int n1 = n % N1;
		const int n2 = n % N2;
		Y[N2 * n1 + n2] = X[n];
	}
	for (int n1 = 0; n1 < N1; n1++) {
		fft_real(Y.data() + n1 * N2, N2);
	}
	for (int n1 = 0; n1 < N1; n1++) {
		X[n1] = Y[N2 * n1];
	}
	fft_real(X, N1);
	for (int n1 = 0; n1 < N1; n1++) {
		Y[N2 * n1] = X[n1];
	}
	for (int k2 = 1; k2 <= N2 / 2; k2++) {
		for (int n1 = 0; n1 < N1; n1++) {
			z[n1].real() = Y[N2 * n1 + k2];
			z[n1].imag() = Y[N2 * n1 - k2 + N2];
		}
		fft(z.data(), N1);
		for (int n1 = 0; n1 < N1; n1++) {
			Y[N2 * n1 + k2] = z[n1].real();
			Y[N2 * n1 - k2 + N2] = z[n1].imag();
		}
	}
	for (int n1 = 0; n1 < N1; n1++) {
		const int n = (N2 * n1) % N;
		X[n] = Y[N2 * n1];
		for (int n2 = 1; n2 <= N2 / 2; n2++) {
			const int n = (N1 * n2 + N2 * n1) % N;
			if (n < N - n) {
				X[n] = Y[N2 * n1 + n2];
				X[N - n] = Y[N2 * n1 - n2 + N2];
			} else {
				X[N - n] = Y[N2 * n1 + n2];
				X[n] = -Y[N2 * n1 - n2 + N2];
			}
		}
	}
	ws.destroy(std::move(Y));
	cws.destroy(std::move(z));
}
