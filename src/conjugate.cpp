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
void fft_conjugate_big(int N1, complex<T>* X, int N) {
	const int N2 = N / N1;
	const int No2 = N / 2;
	const auto& W = twiddles(N);
	const int N1o2 = (N1 + 1) / 2;
	const bool even = N1 % 2 == 0;
	workspace<complex<T>> ws;
	auto z = ws.create(N1);
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
		fft_scramble(z.data(), N1);
		fft(z.data(), N1);
		for (int k1 = 0; k1 < N1; k1++) {
			X[N2 * k1 + k2] = z[k1];
		}
	}
	ws.destroy(std::move(z));
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
	case 5:
		return fft_conjugate<T, 5>(X, N);
	case 6:
		return fft_conjugate<T, 6>(X, N);
	case 7:
		return fft_conjugate<T, 7>(X, N);
	case 8:
		return fft_conjugate<T, 8>(X, N);
	case 9:
		return fft_conjugate<T, 9>(X, N);
	case 10:
		return fft_conjugate<T, 10>(X, N);
	case 11:
		return fft_conjugate<T, 11>(X, N);
	case 12:
		return fft_conjugate<T, 12>(X, N);
	case 13:
		return fft_conjugate<T, 13>(X, N);
	case 14:
		return fft_conjugate<T, 14>(X, N);
	case 15:
		return fft_conjugate<T, 15>(X, N);
	case 16:
		return fft_conjugate<T, 16>(X, N);
	case 17:
		return fft_conjugate<T, 17>(X, N);
	case 18:
		return fft_conjugate<T, 18>(X, N);
	case 19:
		return fft_conjugate<T, 19>(X, N);
	case 20:
		return fft_conjugate<T, 20>(X, N);
	case 21:
		return fft_conjugate<T, 21>(X, N);
	case 22:
		return fft_conjugate<T, 22>(X, N);
	case 23:
		return fft_conjugate<T, 23>(X, N);
	case 24:
		return fft_conjugate<T, 24>(X, N);
	case 25:
		return fft_conjugate<T, 25>(X, N);
	case 26:
		return fft_conjugate<T, 26>(X, N);
	case 27:
		return fft_conjugate<T, 27>(X, N);
	case 28:
		return fft_conjugate<T, 28>(X, N);
	case 29:
		return fft_conjugate<T, 29>(X, N);
	case 30:
		return fft_conjugate<T, 30>(X, N);
	case 31:
		return fft_conjugate<T, 31>(X, N);
	case 32:
		return fft_conjugate<T, 32>(X, N);
	case 33:
		return fft_conjugate<T, 33>(X, N);
	case 34:
		return fft_conjugate<T, 34>(X, N);
	case 35:
		return fft_conjugate<T, 35>(X, N);
	case 36:
		return fft_conjugate<T, 36>(X, N);
	case 37:
		return fft_conjugate<T, 37>(X, N);
	case 38:
		return fft_conjugate<T, 38>(X, N);
	case 39:
		return fft_conjugate<T, 39>(X, N);
	case 40:
		return fft_conjugate<T, 40>(X, N);
	case 41:
		return fft_conjugate<T, 41>(X, N);
	case 42:
		return fft_conjugate<T, 42>(X, N);
	case 43:
		return fft_conjugate<T, 43>(X, N);
	case 44:
		return fft_conjugate<T, 44>(X, N);
	case 45:
		return fft_conjugate<T, 45>(X, N);
	case 46:
		return fft_conjugate<T, 46>(X, N);
	case 47:
		return fft_conjugate<T, 47>(X, N);
	case 48:
		return fft_conjugate<T, 48>(X, N);
	case 49:
		return fft_conjugate<T, 49>(X, N);
	case 50:
		return fft_conjugate<T, 50>(X, N);
	case 51:
		return fft_conjugate<T, 51>(X, N);
	case 52:
		return fft_conjugate<T, 52>(X, N);
	case 53:
		return fft_conjugate<T, 53>(X, N);
	case 54:
		return fft_conjugate<T, 54>(X, N);
	case 55:
		return fft_conjugate<T, 55>(X, N);
	case 56:
		return fft_conjugate<T, 56>(X, N);
	case 57:
		return fft_conjugate<T, 57>(X, N);
	case 58:
		return fft_conjugate<T, 58>(X, N);
	case 59:
		return fft_conjugate<T, 59>(X, N);
	case 60:
		return fft_conjugate<T, 60>(X, N);
	case 61:
		return fft_conjugate<T, 61>(X, N);
	case 62:
		return fft_conjugate<T, 62>(X, N);
	case 63:
		return fft_conjugate<T, 63>(X, N);
	case 64:
		return fft_conjugate<T, 64>(X, N);
	default:
		return fft_conjugate_big<T>(N1, X, N);
	}
}

void fft_conjugate(int N1, complex<double>* X, int N) {
	fft_conjugate1(N1, X, N);
}

void fft_conjugate(int N1, complex<fft_simd4>* X, int N) {
	fft_conjugate1(N1, X, N);
}

