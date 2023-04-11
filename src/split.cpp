#include "fft.hpp"
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
	if (N < FFT_NMAX) {
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



void fft_split(int R, complex<fft_simd4>* X, int N) {
	const int N1 = R;
	const int N1o2 = N1 / 2;
	const int N1o4 = N1 / 4;
	std::vector<complex<fft_simd4>> scratch = create_scratch(2 * N1o2);
	complex<fft_simd4>* ze = scratch.data();
	complex<fft_simd4>* zo = scratch.data() + N1o2;
	const int N2 = N / N1;
	const auto& W = twiddles(N);
	const int No2 = N / 2;
	fft(X, No2);
	for (int n1 = 0; n1 < N1o2; n1++) {
		int o = No2 + N2 * n1;
		fft(X + o, N2);
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
		switch(N1) {
		case 4:
			fft_complex_odd_4((fft_simd4*) zo);
			break;
		case 8:
			fft_complex_odd_8((fft_simd4*) zo);
			break;
		case 16:
			fft_complex_odd_16((fft_simd4*) zo);
			break;
		case 32:
			fft_complex_odd_32((fft_simd4*) zo);
			break;
		}
		for (int k1 = 0; k1 < N1o2; k1++) {
			X[k1 * N2 + k2] = ze[k1] + zo[k1];
			X[k1 * N2 + No2 + k2] = ze[k1] - zo[k1];
		}
	}
	destroy_scratch(std::move(scratch));
}

