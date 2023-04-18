#include "sfft.hpp"
#include "util.hpp"

#include <array>
#include <cmath>
#include <cstring>
#include <numeric>
#include "fftu.hpp"
#include <stack>


void fft_split_conjugate_indices(int R, int* I, int N) {
	if (N < SFFT_NMAX) {
		return;
	}
	const int N1 = R;
	const int N1o2 = N1 / 2;
	const int N1o4 = N1 / 4;
	const int N2 = N / N1;
	const int No2 = N / 2;
	std::vector<int> J;
	J.resize(N);
	for (int n = 0; n < No2; n++) {
		J[n] = I[2 * n];
	}
	for (int n1 = 0; n1 < N1o4; n1++) {
		for (int n2 = 0; n2 < N2; n2++) {
			J[No2 + N2 * n1 + n2] = I[N1 * n2 + 2 * n1 + 1];
		}
	}
	for (int n1 = 0; n1 < N1o4; n1++) {
		for (int n2 = 0; n2 < N2; n2++) {
			J[No2 + N2 * (N1o2 - 1 - n1) + n2] = I[mod(N1 * n2 - 2 * n1 - 1, N)];
		}
	}
	fft_indices(J.data(), No2);
	for (int n1 = 0; n1 < N1o2; n1++) {
		const int o = No2 + N2 * n1;
		fft_indices(J.data() + o, N2);
	}
	std::memcpy(I, J.data(), sizeof(int) * N);
}

template<class T, int N1>
void fft_split_conjugate(complex<T>* X, int N) {
	constexpr int N1o2 = N1 / 2;
	constexpr int N1o4 = N1 / 4;
	std::array<complex<T>, N1o2> ze;
	std::array<complex<T>, N1o2> zo;
	const int N2 = N / N1;
	const auto& W = twiddles(N);
	const int No2 = N / 2;
	fft(X, No2);
	for (int n1 = 0; n1 < N1o2; n1++) {
		int o = No2 + N2 * n1;
		fft(X + o, N2);
	}
	int k2pNo2 = No2;
	for (int k2 = 0; k2 < N2; k2++) {
		for (int n1 = 0; n1 < N1o2; n1++) {
			ze[n1] = X[n1 * N2 + k2];
		}
		int n1m = N1o2 - 1;
		int wi = k2;
		const int twok2 = 2 * k2;
		for (int n1 = 0; n1 < N1o4; n1++) {
			const auto& w = W[wi];
			zo[n1] = X[k2pNo2 + N2 * n1] * w;
			zo[n1m] = X[k2pNo2 + N2 * n1m] * w.conj();
			n1m--;
			wi += twok2;
		}
		sfft_complex_odd<N1>((T*) zo.data());
		int k1N2 = 0;
		for (int k1 = 0; k1 < N1o2; k1++) {
			X[k1N2 + k2] = ze[k1] + zo[k1];
			X[k1N2 + k2pNo2] = ze[k1] - zo[k1];
			k1N2 += N2;
		}
		k2pNo2++;
	}
}

template<class T>
void fft_split_conjugate1(int N1, complex<T>* X, int N) {
	switch (N1) {
	case 4:
		return fft_split_conjugate<T, 4>(X, N);
	case 6:
		return fft_split_conjugate<T, 6>(X, N);
	case 8:
		return fft_split_conjugate<T, 8>(X, N);
	case 10:
		return fft_split_conjugate<T, 10>(X, N);
	case 12:
		return fft_split_conjugate<T, 12>(X, N);
	case 14:
		return fft_split_conjugate<T, 14>(X, N);
	case 16:
		return fft_split_conjugate<T, 16>(X, N);
	case 18:
		return fft_split_conjugate<T, 18>(X, N);
	case 20:
		return fft_split_conjugate<T, 20>(X, N);
	case 22:
		return fft_split_conjugate<T, 22>(X, N);
	case 24:
		return fft_split_conjugate<T, 24>(X, N);
	case 26:
		return fft_split_conjugate<T, 26>(X, N);
	case 28:
		return fft_split_conjugate<T, 28>(X, N);
	case 30:
		return fft_split_conjugate<T, 30>(X, N);
	case 32:
		return fft_split_conjugate<T, 32>(X, N);
	}
}

void fft_split_conjugate(int N1, complex<double>* X, int N) {
	fft_split_conjugate1(N1, X, N);
}

void fft_split_conjugate(int N1, complex<fft_simd4>* X, int N) {
	fft_split_conjugate1(N1, X, N);
}

