#include "fft.hpp"
#include "util.hpp"
#include <cassert>
#include <cstring>

template<int N1, class T>
void fft_cooley_tukey(complex<T>* X, int N) {
	if (N <= SFFT_NMAX) {
		sfft_complex((T*) X, N);
		return;
	}
	const int N2 = N / N1;
	std::array<complex<T>, N1> z;
	const auto& W = twiddles(N);
	for (int n1 = 0; n1 < N1; n1++) {
		fft(X + N2 * n1, N2, false);
	}
	for (int n1 = 0; n1 < N1; n1++) {
		for (int k2 = 0; k2 < N2; k2++) {
			X[n1 * N2 + k2] *= W[n1 * k2];
		}
	}
	for (int k2 = 0; k2 < N2; k2++) {
		for (int n1 = 0; n1 < N1; n1++) {
			z[n1] = X[n1 * N2 + k2];
		}
		sfft_complex((T*) z.data(), N1);
		for (int k1 = 0; k1 < N1; k1++) {
			X[k1 * N2 + k2] = z[k1];
		}
	}
}

template<class T>
void fft_cooley_tukey_select(int N1, complex<T>* X, int N) {
	switch (N1) {
	case 2:
		return fft_cooley_tukey<2, T>(X, N);
	case 3:
		return fft_cooley_tukey<3, T>(X, N);
	case 4:
		return fft_cooley_tukey<4, T>(X, N);
	case 5:
		return fft_cooley_tukey<5, T>(X, N);
	case 6:
		return fft_cooley_tukey<6, T>(X, N);
	case 7:
		return fft_cooley_tukey<7, T>(X, N);
	case 8:
		return fft_cooley_tukey<8, T>(X, N);
	case 9:
		return fft_cooley_tukey<9, T>(X, N);
	case 10:
		return fft_cooley_tukey<10, T>(X, N);
	case 11:
		return fft_cooley_tukey<11, T>(X, N);
	case 12:
		return fft_cooley_tukey<12, T>(X, N);
	case 13:
		return fft_cooley_tukey<13, T>(X, N);
	case 14:
		return fft_cooley_tukey<14, T>(X, N);
	case 15:
		return fft_cooley_tukey<15, T>(X, N);
	case 16:
		return fft_cooley_tukey<16, T>(X, N);
	default:
		assert(false);
		abort();
	}
}

void fft_cooley_tukey(int N1, complex<double>* X, int N) {
	fft_cooley_tukey_select(N1, X, N);
}

void fft_cooley_tukey(int N1, complex<fft_simd4>* X, int N) {
	fft_cooley_tukey_select(N1, X, N);
}

void fft_cooley_tukey_indices(int N1, int* I, int N) {
	if (N <= SFFT_NMAX) {
		return;
	}
	const int N2 = N / N1;
	std::vector<int> J(N);
	for (int n1 = 0; n1 < N1; n1++) {
		for (int n2 = 0; n2 < N2; n2++) {
			J[n1 * N2 + n2] = I[n2 * N1 + n1];
		}
	}
	std::memcpy(I, J.data(), sizeof(int) * N);
	for (int n1 = 0; n1 < N1; n1++) {
		fft_indices(I + N2 * n1, N2);
	}
}

