#include "sfft.hpp"
#include "util.hpp"

#include <array>
#include <cmath>
#include <cstring>
#include <numeric>
#include "fftu.hpp"
#include <stack>

void fft_split_real_indices(int N1, int* I, int N) {
	std::vector<int> J(N);
	const int N2 = N / N1;
	for (int n2 = 0; n2 < N2; n2++) {
		J[n2] = I[N1 * n2 + 0];
		J[n2 + N2] = I[N1 * n2 + 2];
		J[n2 + N / 2] = I[N1 * n2 + 1];
		J[n2 + N / 2 + N2] = I[N1 * n2 + 3];
	}
	std::memcpy(I, J.data(), sizeof(int) * N);
	for (int n1 = 0; n1 < N1; n1++) {
		fft_indices_real(&I[n1 * N2], N2);
	}
}

template<class T, int N1>
void fft_split_real(T* X, int N) {
	if (N <= SFFT_NMAX) {
		sfft_real(X, N);
		return;
	}

	constexpr int N1o2 = N1 / 2;
	const int N2 = N / N1;
	const int No2 = N / 2;
	const int N2p1o2 = (N2 + 1) / 2;

	const auto& W = twiddles(N);

	fft_real(&X[0], N2);
	fft_real(&X[N2], N2);
	fft_real(&X[N / 2], N2);
	fft_real(&X[N / 2 + N2], N2);

	for (int n2 = 0; n2 < N2; n2++) {
		std::swap(X[N2 + n2], X[N / 2 + n2]);
	}

	for (int k2 = 1; k2 < N2p1o2; k2++) {
		for (int n1 = 1; n1 < N1; n1++) {
			complex<T> z;
			const int ir = N2 * n1 + k2;
			const int ii = N2 * n1 - k2 + N2;
			z.real() = X[ir];
			z.imag() = X[ii];
			z *= W[n1 * k2];
			X[ir] = z.real();
			X[ii] = z.imag();
		}
	}

	std::array<T, N1> q;
	for (int n1 = 0; n1 < N1; n1++) {
		q[n1] = X[N2 * n1];
	}
	const auto z0 = q[0] + q[1] + q[2] + q[3];
	const auto z1r = q[0] - q[2];
	const auto z1i = -q[1] + q[3];
	const auto z2 = q[0] - q[1] + q[2] - q[3];
	X[0] = z0;
	X[N2] = z1r;
	X[N / 2] = z2;
	X[N2 + N / 2] = z1i;

	std::array<complex<T>, N1> p;
	for (int k2 = 1; k2 < N2p1o2; k2++) {
		for (int n1 = 0; n1 < N1; n1++) {
			p[n1].real() = X[N2 * n1 + k2];
			p[n1].imag() = X[N2 * n1 - k2 + N2];
		}
		const auto z0 = p[0] + p[1] + p[2] + p[3];
		const auto z1 = p[0] - Ix(p[1]) - p[2] + Ix(p[3]);
		const auto z2 = p[0] - p[1] + p[2] - p[3];
		const auto z3 = p[0] + Ix(p[1]) - p[2] - Ix(p[3]);
		X[N2 * 0 + k2] = z0.real();
		X[N - N2 * 0 - k2] = z0.imag();
		X[N2 * 1 + k2] = z1.real();
		X[N - N2 * 1 - k2] = z1.imag();
		X[N2 * 2 + k2] = -z2.imag();
		X[N - N2 * 2 - k2] = z2.real();
		X[N2 * 3 + k2] = -z3.imag();
		X[N - N2 * 3 - k2] = z3.real();
	}

	if (N2 % 2 == 0) {
		std::array<T, N1> q;
		for (int n1 = 0; n1 < N1; n1++) {
			q[n1] = X[N2 * n1 + N2 / 2];
		}
		const auto t1 = (q[1] - q[3]) * (1.0 / sqrt(2));
		const auto t2 = (q[1] + q[3]) * (1.0 / sqrt(2));
		X[0 + N2 * 0 + N2 / 2] = q[0] + t1;
		X[N - N2 * 0 - N2 / 2] = -q[2] - t2;
		X[0 + N2 * 1 + N2 / 2] = q[0] - t1;
		X[N - N2 * 1 - N2 / 2] = q[2] - t2;
	}
}

template<class T>
void fft_split1_real(int N1, T* X, int N) {
	switch (N1) {
	case 4:
		return fft_split_real<T, 4>(X, N);
	}
}

void fft_split_real(int N1, double* X, int N) {
	fft_split1_real(N1, X, N);
}

void fft_split_real(int N1, fft_simd4* X, int N) {
	fft_split1_real(N1, X, N);
}

