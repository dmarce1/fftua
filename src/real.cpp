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
	const int No2 = N / 2;
	const int N1p1o2 = (N1 + 1) / 2;
	const int N2p1o2 = (N2 + 1) / 2;
	const int N1o2 = N1 / 2;

	const auto& W = twiddles(N);

	std::vector<double> Y;
	Y.resize(N);
	for (int n1 = 0; n1 < N1; n1++) {
		for (int n2 = 0; n2 < N2; n2++) {
			Y[N2 * n1 + n2] = X[N1 * n2 + n1];
		}
	}
	std::memcpy(X, Y.data(), sizeof(double) * N);

	for (int n1 = 0; n1 < N1; n1++) {
		fft_real(N1, &X[n1 * N2], N2);
	}

	for (int k2 = 1; k2 < N2p1o2; k2++) {
		for (int n1 = 1; n1 < N1; n1++) {
			complex<double> z;
			const int ir = N2 * n1 + k2;
			const int ii = N2 * n1 - k2 + N2;
			z.real() = X[ir];
			z.imag() = X[ii];
			z *= W[n1 * k2];
			X[ir] = z.real();
			X[ii] = z.imag();
		}
	}

	{
		double q[N1];
		for (int n1 = 0; n1 < N1; n1++) {
			q[n1] = X[N2 * n1];
		}
		sfft_real(q, N1);
		X[0] = q[0];
		for (int k1 = 1; k1 < N1p1o2; k1++) {
			X[0 + N2 * k1] = q[k1];
			X[N - N2 * k1] = q[N1 - k1];
		}
		if (N1 % 2 == 0) {
			X[No2] = q[N1o2];
		}
	}
	complex<double> p[N1];
	for (int k2 = 1; k2 < N2p1o2; k2++) {
		for (int n1 = 0; n1 < N1; n1++) {
			p[n1].real() = X[N2 * n1 + k2];
			p[n1].imag() = X[N2 * n1 - k2 + N2];
		}
		sfft_complex((double*) p, N1);
		for (int k1 = 0; k1 < N1p1o2; k1++) {
			const int k = N2 * k1 + k2;
			X[k] = p[k1].real();
			X[N - k] = p[k1].imag();
		}
		for (int k1 = N1p1o2; k1 < N1; k1++) {
			const int k = N2 * k1 + k2;
			X[N - k] = p[k1].real();
			X[k] = -p[k1].imag();
		}
	}
	if (N2 % 2 == 0) {
		double q[N1];
		for (int n1 = 0; n1 < N1; n1++) {
			q[n1] = X[N2 * n1 + N2 / 2];
		}
		sfft_skew(q, N1);
		for (int k1 = 0; k1 < N1p1o2; k1++) {
			X[0 + (N2 * k1 + N2 / 2)] = q[k1];
			X[N - (N2 * k1 + N2 / 2)] = q[N1 - k1 - 1];
		}
		if (N1 % 2 == 1) {
			X[No2] = q[N1o2];
		}
	}
}
