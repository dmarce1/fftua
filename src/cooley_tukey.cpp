#include "fftu.hpp"
#include "fft.hpp"
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
		J[n2] = I[N1 * n2];
		J[N / 2 + n2] = I[N1 * n2 + N1 / 2];
		for (int n1 = 1; n1 < N1 / 2; n1++) {
			J[N2 * n1 + n2] = I[N1 * n2 + n1];
		}
		for (int n1 = 1; n1 < N1 / 2; n1++) {
			J[N2 * (N1 - n1) + n2] = I[mod(N1 * n2 - n1, N)];
		}
	}
	for (int n1 = 0; n1 < N1; n1++) {
		fft_indices(J.data() + n1 * N2, N2);
	}
	std::memcpy(I, J.data(), sizeof(int) * N);
}

template<int N1>
void fft_cooley_tukey(complex<fft_simd4>* X, int N) {
	constexpr int N1o2 = N1 / 2;
	const int N2 = N / N1;
	const int No2 = N / 2;
	const auto& W = twiddles(N);
	std::array<complex<fft_simd4>, N1> z;
	for (int n1 = 0; n1 < N1; n1++) {
		fft(X + N2 * n1, N2);
	}
	for (int k2 = 0; k2 < N2; k2++) {
		z[0] = X[k2];
		z[N1o2] = X[No2 + k2] * W[k2 * N1o2];
		auto w = W[k2];
		const auto w0 = w;
		for (int n1 = 1; n1 < N1o2; n1++) {
			z[n1] = X[N2 * n1 + k2] * w;
			z[N1 - n1] = X[N2 * (N1 - n1) + k2] * w.conj();
			w *= w0;
		}
		if (N1 == 2) {
			fft_complex_exe_2((fft_simd4*) z.data());
		} else if (N1 == 4) {
			fft_complex_exe_4((fft_simd4*) z.data());
		} else if (N1 == 8) {
			fft_complex_exe_8((fft_simd4*) z.data());
		} else if (N1 == 16) {
			fft_complex_exe_16((fft_simd4*) z.data());
		} else if (N1 == 32) {
			fft_complex_exe_32((fft_simd4*) z.data());
		}
		for (int k1 = 0; k1 < N1; k1++) {
			X[N2 * k1 + k2] = z[k1];
		}
	}
}

void fft_cooley_tukey(int R, complex<fft_simd4>* X, int N) {
	switch (R) {
	case 2:
		return fft_cooley_tukey<2>(X, N);
	case 4:
		return fft_cooley_tukey<4>(X, N);
	case 8:
		return fft_cooley_tukey<8>(X, N);
	case 16:
		return fft_cooley_tukey<16>(X, N);
	case 32:
		return fft_cooley_tukey<32>(X, N);
	}
}
