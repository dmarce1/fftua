#include "fftu.hpp"
#include "sfft.hpp"
#include "util.hpp"
#include <unordered_map>
#include <cmath>
#include <vector>
#include <stack>
#include <cstring>

void fft_cooley_tukey_indices(int N1, int* I, int N) {
	const int N2 = N / N1;
	std::vector<int> J(N);
	for (int n2 = 0; n2 < N2; n2++) {
		for (int n1 = 0; n1 < N1; n1++) {
			J[N2 * n1 + n2] = I[N1 * n2 + n1];
		}
	}
	for (int n1 = 0; n1 < N1; n1++) {
		fft_indices(J.data() + n1 * N2, N2);
	}
	std::memcpy(I, J.data(), sizeof(int) * N);
}

template<class T, int N1>
void fft_cooley_tukey(complex<T>* X, int N) {
	const int N2 = N / N1;
	const auto& W = twiddles(N);
	std::array<complex<T>, N1> z;
	for (int n1 = 0; n1 < N1; n1++) {
		fft(X + N2 * n1, N2);
	}
	for (int k2 = 0; k2 < N2; k2++) {
		z[0] = X[k2];
		for (int n1 = 1; n1 < N1; n1++) {
			z[n1] = X[N2 * n1 + k2] * W[k2 * n1];
		}
		sfft_complex<N1>((T*) z.data());
		for (int k1 = 0; k1 < N1; k1++) {
			X[N2 * k1 + k2] = z[k1];
		}
	}
}

template<class T>
void fft_cooley_tukey1(int N1, complex<T>* X, int N) {
	switch (N1) {
	case 2:
		return fft_cooley_tukey<T, 2>(X, N);
	case 3:
		return fft_cooley_tukey<T, 3>(X, N);
	case 4:
		return fft_cooley_tukey<T, 4>(X, N);
	case 8:
		return fft_cooley_tukey<T, 8>(X, N);
	case 9:
		return fft_cooley_tukey<T, 9>(X, N);
	case 16:
		return fft_cooley_tukey<T, 16>(X, N);
	case 27:
		return fft_cooley_tukey<T, 27>(X, N);
	case 32:
		return fft_cooley_tukey<T, 32>(X, N);
	case 64:
		return fft_cooley_tukey<T, 64>(X, N);
	}
}

void fft_cooley_tukey(int N1, complex<double>* X, int N) {
	fft_cooley_tukey1(N1, X, N);
}

void fft_cooley_tukey(int N1, complex<fft_simd4>* X, int N) {
	fft_cooley_tukey1(N1, X, N);
}

