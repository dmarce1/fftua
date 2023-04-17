#include "fftu.hpp"
#include "sfft.hpp"
#include "util.hpp"
#include <unordered_map>
#include <cmath>
#include <vector>
#include <stack>
#include <cstring>

void fft_conjugate_indices(int N1, int* I, int N) {
	const int N2 = N / N1;
	std::vector<int> J(N);
	for (int n2 = 0; n2 < N2; n2++) {
		J[n2] = I[N1 * n2];
		if (N1 % 2 == 0) {
			J[N / 2 + n2] = I[N1 * n2 + N1 / 2];
		}
		for (int n1 = 1; n1 < (N1 + 1) / 2; n1++) {
			J[N2 * n1 + n2] = I[N1 * n2 + n1];
		}
		for (int n1 = 1; n1 < (N1 + 1) / 2; n1++) {
			J[N2 * (N1 - n1) + n2] = I[mod(N1 * n2 - n1, N)];
		}
	}
	for (int n1 = 0; n1 < N1; n1++) {
		fft_indices(J.data() + n1 * N2, N2);
	}
	std::memcpy(I, J.data(), sizeof(int) * N);
}

template<class T, int N1>
void fft_conjugate(complex<T>* X, int N) {
	const int N2 = N / N1;
	const int No2 = N / 2;
	const auto& W = twiddles(N);
	constexpr int N1o2 = (N1 + 1) / 2;
	constexpr bool even = N1 % 2 == 0;
	std::array<complex<T>, N1> z;
	for (int n1 = 0; n1 < N1; n1++) {
		fft(X + N2 * n1, N2);
	}
	for (int k2 = 0; k2 < N2; k2++) {
		z[0] = X[k2];
		int n1m = N1 - 1;
		for (int n1 = 1; n1 < N1o2; n1++) {
			const auto& w = W[k2 * n1];
			z[n1] = X[N2 * n1 + k2] * w;
			z[n1m] = X[N2 * n1m + k2] * w.conj();
			n1m--;
		}
		if (even) {
			z[N1o2] = X[No2 + k2] * W[k2 * N1o2];
		}
		sfft_complex<N1>((T*) z.data());
		for (int k1 = 0; k1 < N1; k1++) {
			X[N2 * k1 + k2] = z[k1];
		}
	}
}

template<class T>
void fft_conjugate1(int N1, complex<T>* X, int N) {
	switch (N1) {
	case 2:
		return fft_conjugate<T, 2>(X, N);
	case 3:
		return fft_conjugate<T, 3>(X, N);
	case 4:
		return fft_conjugate<T, 4>(X, N);
	case 8:
		return fft_conjugate<T, 8>(X, N);
	case 9:
		return fft_conjugate<T, 9>(X, N);
	case 16:
		return fft_conjugate<T, 16>(X, N);
	case 27:
		return fft_conjugate<T, 27>(X, N);
	case 32:
		return fft_conjugate<T, 32>(X, N);
	case 64:
		return fft_conjugate<T, 64>(X, N);
	}
}

void fft_conjugate(int N1, complex<double>* X, int N) {
	fft_conjugate1(N1, X, N);
}

void fft_conjugate(int N1, complex<fft_simd4>* X, int N) {
	fft_conjugate1(N1, X, N);
}
