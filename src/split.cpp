#include "fft.hpp"
#include "util.hpp"

#include <array>
#include <cmath>
#include <cstring>
#include <numeric>

static void find_indices(int* I, int N) {
	if (N < FFT_NMAX) {
		return;
	}
	constexpr int N1 = 16;
	constexpr int N1o2 = N1 / 2;
	constexpr int N1o4 = N1 / 4;
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
	find_indices(J.data(), No2);
	for (int n1 = 0; n1 < N1o2; n1++) {
		const int o = No2 + N2 * n1;
		find_indices(J.data() + o, N2);
	}
	std::memcpy(I, J.data(), sizeof(int) * N);
}

template<class T>
static void permute_indices(complex<T>* X, int N) {
	static std::unordered_map<int, std::vector<int>> cache;
	auto iter = cache.find(N);
	if (iter == cache.end()) {
		std::vector<int> I(N);
		std::vector<int> J(N);
		std::iota(I.begin(), I.end(), 0);
		find_indices(I.data(), N);
		for (int n = 0; n < N; n++) {
			J[I[n]] = n;
		}
		cache[N] = std::move(J);
		iter = cache.find(N);
	}
	const auto& I = iter->second;
	static std::vector<bool> flag;
	flag.resize(0);
	flag.resize(N, false);
	for (int n = 0; n < N; n++) {
		if (!flag[n]) {
			flag[n] = true;
			complex<T> tmp = X[n];
			int m = I[n];
			while (m != n) {
				flag[m] = true;
				std::swap(tmp, X[m]);
				m = I[m];
			}
			X[n] = tmp;
		}
	}
}

template<class T>
void fft_split(complex<T>* __restrict__ X, int N, bool permute = true) {
	if (N < FFT_NMAX) {
		fft_complex_simd4((T*) X, N);
		return;
	}
	constexpr int N1 = 16;
	constexpr int N1o2 = N1 / 2;
	constexpr int N1o4 = N1 / 4;
	std::array<complex<T>, N1o2> ze;
	std::array<complex<T>, N1o2> zo;
	const int N2 = N / N1;
	const auto& W = twiddles(N);
	const int No2 = N / 2;
	if (permute) {
		permute_indices(X, N);
	}
	/*	for (int n = 0; n < No2; n++) {
	 Y[n] = X[2 * n];
	 }
	 for (int n1 = 0; n1 < N1o4; n1++) {
	 for (int n2 = 0; n2 < N2; n2++) {
	 Y[No2 + N2 * n1 + n2] = X[N1 * n2 + 2 * n1 + 1];
	 }
	 }
	 for (int n1 = 0; n1 < N1o4; n1++) {
	 for (int n2 = 0; n2 < N2; n2++) {
	 Y[No2 + N2 * (N1o2 - 1 - n1) + n2] = X[mod(N1 * n2 - 2 * n1 - 1, N)];
	 }
	 }*/
	fft_split(X, No2, false);
	for (int n1 = 0; n1 < N1o2; n1++) {
		int o = No2 + N2 * n1;
		fft_split(X + o, N2, false);
	}
	for (int k2 = 0; k2 < N2; k2++) {
		for (int n1 = 0; n1 < N1o2; n1++) {
			ze[n1] = X[n1 * N2 + k2];
		}
		for (int n1 = 0; n1 < N1o4; n1++) {
			const auto& w = W[(2 * n1 + 1) * k2];
			zo[n1] = X[No2 + N2 * n1 + k2] * w;
			zo[N1o2 - 1 - n1] = X[No2 + N2 * (N1o2 - 1 - n1) + k2] * w.conj();
		}
		fft_complex_odd_simd4((T*) zo.data(), N1);
		for (int k1 = 0; k1 < N1o2; k1++) {
			X[k1 * N2 + k2] = ze[k1] + zo[k1];
			X[k1 * N2 + No2 + k2] = ze[k1] - zo[k1];
		}
	}
}

void fft_split_simd(complex<double>* X, complex<double>* Y, int N) {
	if (N < FFT_NMAX) {
		fft_complex((double*) X, N);
		return;
	}
	constexpr int N1 = 4;
	const int No4 = N / 4;
	std::array<complex<fft_simd4>, N1> z;
	const auto& W = vector_twiddles(N1, No4);
	for (int n2 = 0; n2 < No4; n2++) {
		fft_simd4* const y = (fft_simd4*) (Y + n2 * N1);
		const double* const x = (double*) (X + n2 * N1);
		for (int n1 = 0; n1 < N1; n1++) {
			for (int i = 0; i < 2; i++) {
				y[i][n1] = x[2 * n1 + i];
			}
		}
	}
	fft_split((complex<fft_simd4>*) Y,  No4);
	for (int k2 = 0; k2 < No4; k2++) {
		for (int n1 = 0; n1 < N1; n1 += SIMD_SIZE) {
			*((complex<fft_simd4>*) (Y + N1 * k2 + n1)) *= W[n1 / SIMD_SIZE][k2];
		}
	}
	for (int k2 = 0; k2 < No4; k2 += SIMD_SIZE) {
		for (int n1 = 0; n1 < N1; n1++) {
			for (int i = 0; i < SIMD_SIZE; i++) {
				z[n1].real()[i] = ((complex<fft_simd4>*) (Y + N1 * (k2 + i) + (n1 / SIMD_SIZE)))->real()[n1 % SIMD_SIZE];
				z[n1].imag()[i] = ((complex<fft_simd4>*) (Y + N1 * (k2 + i) + (n1 / SIMD_SIZE)))->imag()[n1 % SIMD_SIZE];
			}
		}
		fft_complex_simd4_4((fft_simd4*) z.data());
		for (int k1 = 0; k1 < N1; k1++) {
			for (int i = 0; i < SIMD_SIZE; i++) {
				X[k1 * No4 + k2 + i].real() = z[k1].real()[i];
				X[k1 * No4 + k2 + i].imag() = z[k1].imag()[i];
			}
		}
	}
}
