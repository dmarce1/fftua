#include "fftu.hpp"
#include "util.hpp"
#include <memory>

int reverse_bits(int i, int N) {
	const auto w = ilogb(N);
	int j = 0;
	for (int k = 0; k < w; k++) {
		j <<= 1;
		j |= i & 1;
		i >>= 1;
	}
	return j;
}

int reverse_bits_conj(int i, int N) {
	const auto w = ilogb(N);
	int l = 0;
	int j = i;
	while (l < w) {
		int m = j % 4;
		if (m % 2 == 0) {
			l++;
			j >>= 1;
		} else {
			if (m % 4 == 3) {
				i -= 1 << (l + 2);
			}
			l += 2;
			j >>= 2;
		}
	}
	return reverse_bits(i, N);
}

template<class T>
void scramble(T* x, int N) {
	int N1m = N - 1;
	int j = 0;
	for (int i = 0; i < N1m; i++) {
		if (i > j) {
			std::swap(x[i], x[j]);
		}
		int k = N >> 1;
		while (k <= j) {
			j -= k;
			k >>= 1;
		}
		j += k;
	}
}

#include <numeric>

template<class T>
void fft_width(T* X, int N) {
	const auto& W = twiddles(N);
	constexpr int N1 = 4;
	T T0R, T0I, T1R, T1I, T2R, T2I, T3R, T3I;
	T U0R, U1R, U1I, U2R, U2I, U3R, U3I;
	T C1, C2, C3, S1, S2, S3;
	int i0, i1, i2, i3, i4, i5, i6, i7;
	int j1, j2, j3;
	int N2 = 1;
	while (N2 < N) {
		const int M = N1 * N2;
		const int L = N / M;
		for (int l = 0; l < L; l++) {
			const int L0 = l * M;
			i0 = L0;
			i1 = i0 + N2;
			i2 = i1 + N2;
			i3 = i2 + N2;
			T0R = X[i0] + X[i1];
			T2R = X[i0] - X[i1];
			T1R = X[i2] + X[i3];
			T3R = X[i3] - X[i2];
			X[i0] = T0R + T1R;
			X[i1] = T2R;
			X[i2] = T0R - T1R;
			X[i3] = T3R;
			if (M >= 8) {
				i0 = L0 + N2 / 2;
				i1 = i0 + N2;
				i2 = i1 + N2;
				i3 = i2 + N2;
				U0R = X[i0];
				U2R = X[i1];
				T1R = (X[i2] - X[i3]) * M_SQRT1_2;
				T3R = (X[i2] + X[i3]) * M_SQRT1_2;
				X[i0] = U0R + T1R;
				X[i1] = U0R - T1R;
				X[i2] = U2R - T3R;
				X[i3] = -U2R - T3R;
			}
			for (int k2 = 1; k2 < N2 / 2; k2++) {
				i0 = L0 + k2;
				i1 = i0 + N2;
				i2 = i1 + N2;
				i3 = i2 + N2;
				i4 = L0 + N2 - k2;
				i5 = i4 + N2;
				i6 = i5 + N2;
				i7 = i6 + N2;
				j1 = L * k2;
				j2 = 2 * j1;
				j3 = 3 * j1;
				C1 = W[j1].real();
				C2 = W[j2].real();
				C3 = W[j3].real();
				S1 = W[j1].imag();
				S2 = W[j2].imag();
				S3 = W[j3].imag();
				U1R = X[i2] * C1 - X[i6] * S1;
				U2R = X[i1] * C2 - X[i5] * S2;
				U3R = X[i3] * C3 - X[i7] * S3;
				U1I = X[i2] * S1 + X[i6] * C1;
				U2I = X[i1] * S2 + X[i5] * C2;
				U3I = X[i3] * S3 + X[i7] * C3;
				T0R = X[i0] + U2R;
				T0I = X[i4] + U2I;
				T2R = X[i0] - U2R;
				T2I = X[i4] - U2I;
				T1R = U1R + U3R;
				T1I = U1I + U3I;
				T3R = U3R - U1R;
				T3I = U1I - U3I;
				X[i0] = T0R + T1R;
				X[i1] = T2R + T3I;
				X[i2] = T1I - T0I;
				X[i3] = T3R - T2I;
				X[i4] = T2R - T3I;
				X[i5] = T0R - T1R;
				X[i6] = T2I + T3R;
				X[i7] = T0I + T1I;
			}
		}
		N2 *= N1;
	}
}

template<class T>
void fft_depth(T* X, int N) {
	constexpr int N1 = 4;
	const int N2 = N / N1;

	T T0R, T0I, T1R, T1I, T2R, T2I, T3R, T3I;
	T U0R, U1R, U1I, U2R, U2I, U3R, U3I;
	T C1, C2, C3, S1, S2, S3;
	int i0, i1, i2, i3, i4, i5, i6, i7;
	int j1, j2, j3;

	if (N <= 1) {
		return;
	} else if (N == 2) {
		U0R = X[0];
		U1R = X[1];
		X[0] = U0R + U1R;
		X[1] = U0R - U1R;
		return;
	}

	for (int n1 = 0; n1 < N1; n1++) {
		fft_depth(&X[n1 * N2], N2);
	}

	i0 = 0;
	i1 = i0 + N2;
	i2 = i1 + N2;
	i3 = i2 + N2;
	T0R = X[i0] + X[i1];
	T2R = X[i0] - X[i1];
	T1R = X[i2] + X[i3];
	T3R = X[i3] - X[i2];
	X[i0] = T0R + T1R;
	X[i1] = T2R;
	X[i2] = T0R - T1R;
	X[i3] = T3R;
	if (N >= 8) {
		i0 = N2 / 2;
		i1 = i0 + N2;
		i2 = i1 + N2;
		i3 = i2 + N2;
		U0R = X[i0];
		U2R = X[i1];
		T1R = (X[i2] - X[i3]) * M_SQRT1_2;
		T3R = (X[i2] + X[i3]) * M_SQRT1_2;
		X[i0] = U0R + T1R;
		X[i1] = U0R - T1R;
		X[i2] = U2R - T3R;
		X[i3] = -U2R - T3R;
	}
	const auto& W = twiddles(N);
	for (int k2 = 1; k2 < N2 / 2; k2++) {
		i0 = k2;
		i1 = i0 + N2;
		i2 = i1 + N2;
		i3 = i2 + N2;
		i4 = N2 - k2;
		i5 = i4 + N2;
		i6 = i5 + N2;
		i7 = i6 + N2;
		j1 = k2;
		j2 = 2 * j1;
		j3 = 3 * j1;
		C1 = W[j1].real();
		C2 = W[j2].real();
		C3 = W[j3].real();
		S1 = W[j1].imag();
		S2 = W[j2].imag();
		S3 = W[j3].imag();
		U1R = X[i2] * C1 - X[i6] * S1;
		U2R = X[i1] * C2 - X[i5] * S2;
		U3R = X[i3] * C3 - X[i7] * S3;
		U1I = X[i2] * S1 + X[i6] * C1;
		U2I = X[i1] * S2 + X[i5] * C2;
		U3I = X[i3] * S3 + X[i7] * C3;
		T0R = X[i0] + U2R;
		T0I = X[i4] + U2I;
		T2R = X[i0] - U2R;
		T2I = X[i4] - U2I;
		T1R = U1R + U3R;
		T1I = U1I + U3I;
		T3R = U3R - U1R;
		T3I = U1I - U3I;
		X[i0] = T0R + T1R;
		X[i1] = T2R + T3I;
		X[i2] = T1I - T0I;
		X[i3] = T3R - T2I;
		X[i4] = T2R - T3I;
		X[i5] = T0R - T1R;
		X[i6] = T2I + T3R;
		X[i7] = T0I + T1I;
	}
}

void fft_2pow(double* x, int N) {
	scramble((fft_simd4*) x, N / 4);
	fft_depth((fft_simd4*) x, N / 4);
}

