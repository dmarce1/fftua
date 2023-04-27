#include "sfft.hpp"
#include "util.hpp"

#include <array>
#include <cmath>
#include <cstring>
#include <numeric>
#include "fftu.hpp"
#include <stack>

int bsqrt(int N) {
	auto l = ilogb(N);
	l /= 2;
	return 1 << l;
}

void fft_split_indices(int R, int* I, int N) {
	if (N <= SFFT_NMAX) {
		return;
	}
	const int N1 = R;
	const int N1o2 = N1 / 2;
	const int N2 = N / N1;
	const int No2 = N / 2;
	std::vector<int> J;
	J.resize(N);
	for (int n = 0; n < No2; n++) {
		J[n] = I[2 * n];
	}
	for (int n1 = 0; n1 < N1o2; n1++) {
		for (int n2 = 0; n2 < N2; n2++) {
			J[No2 + N2 * n1 + n2] = I[N1 * n2 + 2 * n1 + 1];
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
void fft_split(complex<T>* X, int N) {
	constexpr int N1o2 = N1 / 2;
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
		int wi = k2;
		const int twok2 = 2 * k2;
		for (int n1 = 0; n1 < N1o2; n1++) {
			const auto& w = W[wi];
			zo[n1] = X[k2pNo2 + N2 * n1] * w;
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
void fft_split1(int N1, complex<T>* X, int N) {
	switch (N1) {
	case 4:
		return fft_split<T, 4>(X, N);
	case 8:
		return fft_split<T, 8>(X, N);
	case 12:
		return fft_split<T, 12>(X, N);
	case 16:
		return fft_split<T, 16>(X, N);
	case 20:
		return fft_split<T, 20>(X, N);
	case 24:
		return fft_split<T, 24>(X, N);
	case 28:
		return fft_split<T, 28>(X, N);
	case 32:
		return fft_split<T, 32>(X, N);
	case 36:
		return fft_split<T, 36>(X, N);
	case 40:
		return fft_split<T, 40>(X, N);
	case 44:
		return fft_split<T, 44>(X, N);
	case 48:
		return fft_split<T, 48>(X, N);
	case 52:
		return fft_split<T, 52>(X, N);
	case 56:
		return fft_split<T, 56>(X, N);
	case 60:
		return fft_split<T, 60>(X, N);
	case 64:
		return fft_split<T, 64>(X, N);
	}
}

void fft_split(int N1, complex<double>* X, int N) {
	fft_split1(N1, X, N);
}

void fft_split(int N1, complex<fft_simd4>* X, int N) {
	fft_split1(N1, X, N);
}

