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
void fft_cooley_tukey_big(int N1, complex<T>* X, int N) {
	const int N2 = N / N1;
	//printf( "%i %i\n", N1, N2);
	const auto& W = twiddles(N);
	workspace<complex<T>> ws;
	auto z = ws.create(N1);
	for (int n1 = 0; n1 < N1; n1++) {
		fft(X + N2 * n1, N2);
	}
	for (int k2 = 0; k2 < N2; k2++) {
		z[0] = X[k2];
		for (int n1 = 1; n1 < N1; n1++) {
			z[n1] = X[N2 * n1 + k2] * W[k2 * n1];
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
void fft_cooley_tukey1(int N1, complex<T>* X, int N) {
	switch (N1) {
	case 2:
		return fft_cooley_tukey<T, 2>(X, N);
	case 3:
		return fft_cooley_tukey<T, 3>(X, N);
	case 4:
		return fft_cooley_tukey<T, 4>(X, N);
	case 5:
		return fft_cooley_tukey<T, 5>(X, N);
	case 6:
		return fft_cooley_tukey<T, 6>(X, N);
	case 7:
		return fft_cooley_tukey<T, 7>(X, N);
	case 8:
		return fft_cooley_tukey<T, 8>(X, N);
	case 9:
		return fft_cooley_tukey<T, 9>(X, N);
	case 10:
		return fft_cooley_tukey<T, 10>(X, N);
	case 11:
		return fft_cooley_tukey<T, 11>(X, N);
	case 12:
		return fft_cooley_tukey<T, 12>(X, N);
	case 13:
		return fft_cooley_tukey<T, 13>(X, N);
	case 14:
		return fft_cooley_tukey<T, 14>(X, N);
	case 15:
		return fft_cooley_tukey<T, 15>(X, N);
	case 16:
		return fft_cooley_tukey<T, 16>(X, N);
	case 17:
		return fft_cooley_tukey<T, 17>(X, N);
	case 18:
		return fft_cooley_tukey<T, 18>(X, N);
	case 19:
		return fft_cooley_tukey<T, 19>(X, N);
	case 20:
		return fft_cooley_tukey<T, 20>(X, N);
	case 21:
		return fft_cooley_tukey<T, 21>(X, N);
	case 22:
		return fft_cooley_tukey<T, 22>(X, N);
	case 23:
		return fft_cooley_tukey<T, 23>(X, N);
	case 24:
		return fft_cooley_tukey<T, 24>(X, N);
	case 25:
		return fft_cooley_tukey<T, 25>(X, N);
	case 26:
		return fft_cooley_tukey<T, 26>(X, N);
	case 27:
		return fft_cooley_tukey<T, 27>(X, N);
	case 28:
		return fft_cooley_tukey<T, 28>(X, N);
	case 29:
		return fft_cooley_tukey<T, 29>(X, N);
	case 30:
		return fft_cooley_tukey<T, 30>(X, N);
	case 31:
		return fft_cooley_tukey<T, 31>(X, N);
	case 32:
		return fft_cooley_tukey<T, 32>(X, N);
	case 33:
		return fft_cooley_tukey<T, 33>(X, N);
	case 34:
		return fft_cooley_tukey<T, 34>(X, N);
	case 35:
		return fft_cooley_tukey<T, 35>(X, N);
	case 36:
		return fft_cooley_tukey<T, 36>(X, N);
	case 37:
		return fft_cooley_tukey<T, 37>(X, N);
	case 38:
		return fft_cooley_tukey<T, 38>(X, N);
	case 39:
		return fft_cooley_tukey<T, 39>(X, N);
	case 40:
		return fft_cooley_tukey<T, 40>(X, N);
	case 41:
		return fft_cooley_tukey<T, 41>(X, N);
	case 42:
		return fft_cooley_tukey<T, 42>(X, N);
	case 43:
		return fft_cooley_tukey<T, 43>(X, N);
	case 44:
		return fft_cooley_tukey<T, 44>(X, N);
	case 45:
		return fft_cooley_tukey<T, 45>(X, N);
	case 46:
		return fft_cooley_tukey<T, 46>(X, N);
	case 47:
		return fft_cooley_tukey<T, 47>(X, N);
	case 48:
		return fft_cooley_tukey<T, 48>(X, N);
	case 49:
		return fft_cooley_tukey<T, 49>(X, N);
	case 50:
		return fft_cooley_tukey<T, 50>(X, N);
	case 51:
		return fft_cooley_tukey<T, 51>(X, N);
	case 52:
		return fft_cooley_tukey<T, 52>(X, N);
	case 53:
		return fft_cooley_tukey<T, 53>(X, N);
	case 54:
		return fft_cooley_tukey<T, 54>(X, N);
	case 55:
		return fft_cooley_tukey<T, 55>(X, N);
	case 56:
		return fft_cooley_tukey<T, 56>(X, N);
	case 57:
		return fft_cooley_tukey<T, 57>(X, N);
	case 58:
		return fft_cooley_tukey<T, 58>(X, N);
	case 59:
		return fft_cooley_tukey<T, 59>(X, N);
	case 60:
		return fft_cooley_tukey<T, 60>(X, N);
	case 61:
		return fft_cooley_tukey<T, 61>(X, N);
	case 62:
		return fft_cooley_tukey<T, 62>(X, N);
	case 63:
		return fft_cooley_tukey<T, 63>(X, N);
	case 64:
		return fft_cooley_tukey<T, 64>(X, N);
	default:
		return fft_cooley_tukey_big<T>(N1, X, N);
	}
}

void fft_cooley_tukey(int N1, complex<double>* X, int N) {
	fft_cooley_tukey1(N1, X, N);
}

void fft_cooley_tukey(int N1, complex<fft_simd4>* X, int N) {
	fft_cooley_tukey1(N1, X, N);
}

