#include "sfft.hpp"
#include "util.hpp"

#include <array>
#include <cmath>
#include <cstring>
#include <numeric>
#include "fftu.hpp"
#include <stack>

void fft_split_real_indices(int N1, int* I, int N) {
}

template<class T, int N1>
void fft_split_real(T* X, int N) {
	if (N <= SFFT_NMAX) {
		sfft_real(X, N);
		return;
	}
	X--;
	int j = 1;
	for (int i = 1; i < N; i++) {
		if (i < j) {
			std::swap(X[i], X[j]);
		}
		int k = N / 2;
		while (k < j) {
			j -= k;
			k /= 2;
		}
		j += k;
	}
	int is = 1;
	int id = 4;
	do {
		for (int i0 = is; i0 <= N; i0 += id) {
			int i1 = i0 + 1;
			T R1 = X[i0];
			X[i0] = R1 + X[i1];
			X[i1] = R1 - X[i1];
		}
		is = 2 * id - 1;
		id = 4 * id;
	} while (is < N);
	int N2 = 2;
	int M = ilogb(N);
	for (int k = 2; k <= M; k++) {
		N2 *= 2;
		int N4 = N2 / 4;
		int N8 = N2 / 8;
		const double E = 2.0 * M_PI / N2;
		is = 0;
		id = N2 * 2;
		while (is < N) {
			for (int i = is; i < N; i += id) {
				int i1 = i + 1;
				int i2 = i1 + N4;
				int i3 = i2 + N4;
				int i4 = i3 + N4;
				T T1 = X[i4] + X[i3];
				T T2;
				X[i4] = X[i4] - X[i3];
				X[i3] = X[i1] - T1;
				X[i1] = X[i1] + T1;
				if (N4 != 1) {
					i1 += N8;
					i2 += N8;
					i3 += N8;
					i4 += N8;
					T1 = (X[i3] + X[i4]) * (1.0 / sqrt(2.0));
					T2 = (X[i3] - X[i4]) * (1.0 / sqrt(2.0));
					X[i4] = X[i2] - T1;
					X[i3] = -X[i2] - T1;
					X[i2] = X[i1] - T2;
					X[i1] = X[i1] + T2;
				}
			}
			is = 2 * id - N2;
			id = 4 * id;
		}
		double A = E;
		for (int j = 2; j <= N8; j++) {
			double A3 = 3.0 * A;
			double CC1 = cos(A);
			double SS1 = sin(A);
			double CC3 = cos(A3);
			double SS3 = sin(A3);
			A = j * E;
			is = 0;
			id = 2 * N2;
			while (is < N) {
				for (int i = is; i < N; i += id) {
					int i1 = i + j;
					int i2 = i1 + N4;
					int i3 = i2 + N4;
					int i4 = i3 + N4;
					int i5 = i + N4 - j + 2;
					int i6 = i5 + N4;
					int i7 = i6 + N4;
					int i8 = i7 + N4;
					T T1 = X[i3] * CC1 + X[i7] * SS1;
					T T2 = X[i7] * CC1 - X[i3] * SS1;
					T T3 = X[i4] * CC3 + X[i8] * SS3;
					T T4 = X[i8] * CC3 - X[i4] * SS3;
					T T5 = T1 + T3;
					T T6 = T2 + T4;
					T3 = T1 - T3;
					T4 = T2 - T4;
					T2 = X[i6] + T6;
					X[i3] = T6 - X[i6];
					X[i8] = T2;
					T2 = X[i2] - T3;
					X[i7] = -X[i2] - T3;
					X[i4] = T2;
					T1 = X[i1] + T5;
					X[i6] = X[i1] - T5;
					X[i1] = T1;
					T1 = X[i5] + T4;
					X[i5] = X[i5] - T4;
					X[i2] = T1;
				}
				is = 2 * id - N2;
				id = 4 * id;
			}
		}
	}
}

template<class T>
void fft_split1_real(int N1, T* X, int N) {
	switch (N1) {
	case 4:
		return fft_split_real<T, 4>(X, N);
	}
}

void fft_split_real(int N1, double* X, int N) {
	fft_split1_real(N1, X, N);
}

void fft_split_real(int N1, fft_simd4* X, int N) {
	fft_split1_real(N1, X, N);
}

