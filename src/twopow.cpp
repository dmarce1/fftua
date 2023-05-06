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

template<class T>
void scramble_into(T* x, T* y, int N) {
	int N1m = N - 1;
	int j = 0;
	for (int i = 0; i < N1m; i++) {
		y[j] = x[i];
		int k = N >> 1;
		while (k <= j) {
			j -= k;
			k >>= 1;
		}
		j += k;
	}
	y[N1m] = x[N1m];
}

#include <numeric>

#include <numeric>

template<class T>
void fft_depth(T* X, int N) {
	constexpr int N1 = 4;
	if (N <= SFFT_NMAX) {
		scramble(X, N);
		sfft_real(X, N);
		return;
	}

	constexpr int N1o2 = N1 / 2;
	const int N2 = N / N1;
	const int No2 = N / 2;
	const int N2o2 = N2 / 2;

	static std::vector<int> D;
	if (!D.size()) {
		D.resize(N1);
		std::iota(D.begin(), D.end(), 0);
		scramble(D.data(), N1);
	}

	const auto& W = twiddles(N);

	for (int n1 = 0; n1 < N1; n1++) {
		fft_depth(&X[n1 * N2], N2);
	}

	std::array<complex<T>, N1> p;
	std::array<T, N1> q;
	for (int n1 = 0; n1 < N1; n1++) {
		q[n1] = X[N2 * D[n1]];
	}
	sfft_real(q.data(), N1);
	X[0] = q[0];
	X[No2] = q[N1o2];
	for (int k1 = 1; k1 < N1o2; k1++) {
		X[N2 * k1] = q[k1];
		X[N - N2 * k1] = q[N1 - k1];
	}
	if (N2o2) {
		for (int n1 = 0; n1 < N1; n1++) {
			q[n1] = X[N2 * D[n1] + N2o2];
		}
		sfft_skew(q.data(), N1);
		for (int k1 = 0; k1 < N1o2; k1++) {
			X[(N2 * k1 + N2 / 2)] = q[k1];
			X[N - (N2 * k1 + N2 / 2)] = q[N1 - k1 - 1];
		}
	}
	for (int k2 = 1; k2 < N2o2; k2++) {
		for (int n1 = 0; n1 < N1; n1++) {
			p[n1].real() = X[N2 * D[n1] + k2];
			p[n1].imag() = X[N2 * D[n1] - k2 + N2];
		}
		for (int n1 = 1; n1 < N1; n1++) {
			p[n1] *= W[n1 * k2];
		}
		sfft_complex((T*) p.data(), N1);
		for (int k1 = 0; k1 < N1o2; k1++) {
			const int k = N2 * k1 + k2;
			X[k] = p[k1].real();
			X[N - k] = p[k1].imag();
		}
		for (int k1 = N1o2; k1 < N1; k1++) {
			const int k = N2 * k1 + k2;
			X[N - k] = p[k1].real();
			X[k] = -p[k1].imag();
		}
	}
}

template<class T>
void fft_complex(T* X, T* Y, int N) {
	constexpr int N1 = 4;
	T T0R, T0I, T1R, T1I, T2R, T2I, T3R, T3I;
	T U0R, U0I, U1R, U1I, U2R, U2I, U3R, U3I;
	T C1, C2, C3, S1, S2, S3;
	int i0, i1, i2, i3;
	int j1, j2, j3;
	const int D = N / N1;
	int N2 = N / N1;
	int M1 = ilogb(N) - ilogb(N1);
	int M2 = M1 - ilogb(N1);
	const auto& W = twiddles(N);
	while (N2 >= 1) {
		const int M = N1 * N2;
		const int L = N / M;
		for (int k2 = 0; k2 < N / N1; k2++) {
			i0 = k2;
			i1 = i0 + D;
			i2 = i1 + D;
			i3 = i2 + D;
			U0R = X[i0];
			U1R = X[i1];
			U2R = X[i2];
			U3R = X[i3];
			U0I = Y[i0];
			U1I = Y[i1];
			U2I = Y[i2];
			U3I = Y[i3];
			T0R = U0R + U2R;
			T0I = U0I + U2I;
			T2R = U0R - U2R;
			T2I = U0I - U2I;
			T1R = U1R + U3R;
			T1I = U1I + U3I;
			T3R = U1R - U3R;
			T3I = U1I - U3I;
			U0R = T0R + T1R;
			U0I = T0I + T1I;
			U1R = T2R + T3I;
			U1I = T2I - T3R;
			U2R = T0R - T1R;
			U2I = T0I - T1I;
			U3R = T2R - T3I;
			U3I = T2I + T3R;
			j1 = L * (k2 % (M / N1));
			j2 = 2 * j1;
			j3 = 3 * j1;
			C1 = W[j1].real();
			C2 = W[j2].real();
			C3 = W[j3].real();
			S1 = W[j1].imag();
			S2 = W[j2].imag();
			S3 = W[j3].imag();
			T1R = U1R;
			T2R = U2R;
			T3R = U3R;
			U1R = T1R * C1 - U1I * S1;
			U2R = T2R * C2 - U2I * S2;
			U3R = T3R * C3 - U3I * S3;
			U1I = T1R * S1 + U1I * C1;
			U2I = T2R * S2 + U2I * C2;
			U3I = T3R * S3 + U3I * C3;
			X[i0] = U0R;
			Y[i0] = U0I;
			X[i1] = U1R;
			Y[i1] = U1I;
			X[i2] = U2R;
			Y[i2] = U2I;
			X[i3] = U3R;
			Y[i3] = U3I;
		}
		const int mask1 = (1 << M1) | (2 << M1);
		const int mask2 = (1 << M2) | (2 << M2);
		if (M2 >= 0) {
			for (int i = 0; i < N; i++) {
				int j = i;
				int p1 = (j & mask1) >> M1;
				int p2 = (j & mask2) >> M2;
				j &= ~mask1;
				j &= ~mask2;
				j |= p1 << M2;
				j |= p2 << M1;
				if (i < j) {
					std::swap(X[i], X[j]);
					std::swap(Y[i], Y[j]);
				}
			}
		}
		M2 -= ilogb(N1);
		N2 /= N1;
	}
	for (int n = 0; n < N1; n++) {
		int N2 = N / N1;
		auto* x = X + n * N / N1;
		auto* y = Y + n * N / N1;
		int N1m = N2 - 1;
		int j = 0;
		for (int i = 0; i < N1m; i++) {
			if (i > j) {
				std::swap(x[i], x[j]);
				std::swap(y[i], y[j]);
			}
			int k = N2 / 4;
			while (3 * k <= j) {
				j -= 3 * k;
				k >>= 2;
			}
			j += k;
		}
	}
}

template<class T>
void fft_complex2(T* X, T* Y, int N) {
	constexpr int N1 = 2;
	T T0R, T0I, T1R, T1I;
	T U0R, U0I, U1R, U1I;
	T C1, S1;
	int i0, i1;
	int j1;
	const int D = N / N1;
	int N2 = N / N1;
	int M1 = ilogb(N) - ilogb(N1);
	int M2 = M1 - ilogb(N1);
	while (N2 >= 1) {
		const int M = N1 * N2;
		for (int k2 = 0; k2 < N / N1; k2++) {
			const auto& W = twiddles(M);
			i0 = k2;
			i1 = i0 + D;
			U0R = X[i0];
			U1R = X[i1];
			U0I = Y[i0];
			U1I = Y[i1];
			T0R = U0R + U1R;
			T0I = U0I + U1I;
			T1R = U0R - U1R;
			T1I = U0I - U1I;
			U0R = T0R;
			U0I = T0I;
			U1R = T1R;
			U1I = T1I;
			j1 = k2 % (M / N1);
			C1 = W[j1].real();
			S1 = W[j1].imag();
			T1R = U1R;
			U1R = T1R * C1 - U1I * S1;
			U1I = T1R * S1 + U1I * C1;
			X[i0] = U0R;
			Y[i0] = U0I;
			X[i1] = U1R;
			Y[i1] = U1I;
		}
		const int mask1 = (1 << M1);
		const int mask2 = (1 << M2);
		if (M2 >= 0) {
			for (int i = 0; i < N; i++) {
				int j = i;
				int p1 = (j & mask1) >> M1;
				int p2 = (j & mask2) >> M2;
				j &= ~mask1;
				j &= ~mask2;
				j |= p1 << M2;
				j |= p2 << M1;
				if (i < j) {
					std::swap(X[i], X[j]);
					std::swap(Y[i], Y[j]);
				}
			}
		}
		M2 -= ilogb(N1);
		N2 /= N1;
	}
	for (int n = 0; n < 2; n++) {
		scramble(X + n * N / 2, N / 2);
		scramble(Y + n * N / 2, N / 2);
	}
}

template<class T>
void fft_dit(T* x, T* y, int N) {
	constexpr int N1 = 8;
	const int N2 = N / N1;
	T T0R, T0I, T1R, T1I, T2R, T2I, T3R, T3I;
	T T4R, T4I, T5R, T5I, T6R, T6I, T7R, T7I;
	T U0R, U0I, U1R, U1I, U2R, U2I, U3R, U3I;
	T U4R, U4I, U5R, U5I, U6R, U6I, U7R, U7I;
	T C1, C2, C3, C4, C5, C6, C7, S1, S2, S3, S4, S5, S6, S7;
	int i0, i1, i2, i3, i4, i5, i6, i7;
	int j1, j2, j3, j4, j5, j6, j7;

	if (N <= 1) {
		return;
	} else if (N == 2) {
		T0R = x[0];
		T0I = y[0];
		T1R = x[1];
		T1I = y[1];
		U0R = T0R + T1R;
		U0I = T0I + T1I;
		U1R = T0R - T1R;
		U1I = T0I - T1I;
		x[0] = U0R;
		x[1] = U1R;
		y[0] = U0I;
		y[1] = U1I;
		return;
	} else if (N == 4) {
		U0R = x[0];
		U1R = x[2];
		U2R = x[1];
		U3R = x[3];
		U0I = y[0];
		U1I = y[2];
		U2I = y[1];
		U3I = y[3];
		T0R = U0R + U2R;
		T0I = U0I + U2I;
		T2R = U0R - U2R;
		T2I = U0I - U2I;
		T1R = U1R + U3R;
		T1I = U1I + U3I;
		T3R = U1R - U3R;
		T3I = U1I - U3I;
		U0R = T0R + T1R;
		U0I = T0I + T1I;
		U1R = T2R + T3I;
		U1I = T2I - T3R;
		U2R = T0R - T1R;
		U2I = T0I - T1I;
		U3R = T2R - T3I;
		U3I = T2I + T3R;
		x[0] = U0R;
		y[0] = U0I;
		x[1] = U1R;
		y[1] = U1I;
		x[2] = U2R;
		y[2] = U2I;
		x[3] = U3R;
		y[3] = U3I;
	}

	for (int n1 = 0; n1 < N1; n1++) {
		fft_dit(x + n1 * N2, y + n1 * N2, N2);
	}
	const auto& W = twiddles(N);
	for (int k2 = 0; k2 < N2; k2++) {
		j1 = k2;
		j2 = 2 * k2;
		j3 = 3 * k2;
		j4 = 4 * k2;
		j5 = 5 * k2;
		j6 = 6 * k2;
		j7 = 7 * k2;
		C1 = W[j1].real();
		C2 = W[j2].real();
		C3 = W[j3].real();
		C4 = W[j4].real();
		C5 = W[j5].real();
		C6 = W[j6].real();
		C7 = W[j7].real();
		S1 = W[j1].imag();
		S2 = W[j2].imag();
		S3 = W[j3].imag();
		S4 = W[j4].imag();
		S5 = W[j5].imag();
		S6 = W[j6].imag();
		S7 = W[j7].imag();
		i0 = k2;
		i1 = i0 + N2;
		i2 = i1 + N2;
		i3 = i2 + N2;
		i4 = i3 + N2;
		i5 = i4 + N2;
		i6 = i5 + N2;
		i7 = i6 + N2;
		U0R = x[i0];
		U1R = x[i4];
		U2R = x[i2];
		U3R = x[i6];
		U4R = x[i1];
		U5R = x[i5];
		U6R = x[i3];
		U7R = x[i7];
		U0I = y[i0];
		U1I = y[i4];
		U2I = y[i2];
		U3I = y[i6];
		U4I = y[i1];
		U5I = y[i5];
		U6I = y[i3];
		U7I = y[i7];
		T1R = U1R;
		T2R = U2R;
		T3R = U3R;
		T4R = U4R;
		T5R = U5R;
		T6R = U6R;
		T7R = U7R;
		U1R = T1R * C1 - U1I * S1;
		U2R = T2R * C2 - U2I * S2;
		U3R = T3R * C3 - U3I * S3;
		U4R = T4R * C4 - U4I * S4;
		U5R = T5R * C5 - U5I * S5;
		U6R = T6R * C6 - U6I * S6;
		U7R = T7R * C7 - U7I * S7;
		U1I = T1R * S1 + U1I * C1;
		U2I = T2R * S2 + U2I * C2;
		U3I = T3R * S3 + U3I * C3;
		U4I = T4R * S4 + U4I * C4;
		U5I = T5R * S5 + U5I * C5;
		U6I = T6R * S6 + U6I * C6;
		U7I = T7R * S7 + U7I * C7;
		T0R = U0R + U4R;
		T0I = U0I + U4I;
		T2R = U0R - U4R;
		T2I = U0I - U4I;
		T1R = U2R + U6R;
		T1I = U2I + U6I;
		T3R = U2R - U6R;
		T3I = U2I - U6I;
		T4R = U1R + U5R;
		T4I = U1I + U5I;
		T6R = U1R - U5R;
		T6I = U1I - U5I;
		T5R = U3R + U7R;
		T5I = U3I + U7I;
		T7R = U3R - U7R;
		T7I = U3I - U7I;
		U0R = T0R + T1R;
		U0I = T0I + T1I;
		U2R = T2R + T3I;
		U2I = T2I - T3R;
		U4R = T0R - T1R;
		U4I = T0I - T1I;
		U6R = T2R - T3I;
		U6I = T2I + T3R;
		U1R = T4R + T5R;
		U1I = T4I + T5I;
		U3R = T6R + T7I;
		U3I = T6I - T7R;
		U5R = T4R - T5R;
		U5I = T4I - T5I;
		U7R = T6R - T7I;
		U7I = T6I + T7R;
		T3R = U3R;
		T5R = U5R;
		T7R = U7R;
		U3R = (T3R + U3I) * M_SQRT1_2;
		U3I = (-T3R + U3I) * M_SQRT1_2;
		U5R = U5I;
		U5I = -T5R;
		U7R = (-T7R + U7I) * M_SQRT1_2;
		U7I = (-T7R - U7I) * M_SQRT1_2;
		T0R = U0R;
		T0I = U0I;
		T2R = U2R;
		T2I = U2I;
		T4R = U4R;
		T4I = U4I;
		T6R = U6R;
		T6I = U6I;
		U0R = T0R + U1R;
		U0I = T0I + U1I;
		U1R = T0R - U1R;
		U1I = T0I - U1I;
		U2R = T2R + U3R;
		U2I = T2I + U3I;
		U3R = T2R - U3R;
		U3I = T2I - U3I;
		U4R = T4R + U5R;
		U4I = T4I + U5I;
		U5R = T4R - U5R;
		U5I = T4I - U5I;
		U6R = T6R + U7R;
		U6I = T6I + U7I;
		U7R = T6R - U7R;
		U7I = T6I - U7I;
		x[i0] = U0R;
		y[i0] = U0I;
		x[i1] = U2R;
		y[i1] = U2I;
		x[i2] = U4R;
		y[i2] = U4I;
		x[i3] = U6R;
		y[i3] = U6I;
		x[i4] = U1R;
		y[i4] = U1I;
		x[i5] = U3R;
		y[i5] = U3I;
		x[i6] = U5R;
		y[i6] = U5I;
		x[i7] = U7R;
		y[i7] = U7I;
	}
}

void fft_2pow_complex(double* x, int N) {
	std::vector<double> re(N), im(N);
	for (int n = 0; n < N; n++) {
		re[n] = x[2 * n];
		im[n] = x[2 * n + 1];
	}

	scramble(re.data(), N);
	scramble(im.data(), N);
	fft_dit(re.data(), im.data(), N);

	for (int n = 0; n < N; n++) {
		x[2 * n] = re[n];
		x[2 * n + 1] = im[n];
	}
}

void fft_2pow(double* x, int N) {
	static std::vector<double> y;
	constexpr int N1 = 4;
	const int N2 = N / N1;
	y.resize(N);
	scramble_into((fft_simd4*) x, (fft_simd4*) y.data(), N2);
	fft_depth((fft_simd4*) y.data(), N2);
	const auto& W = vector_twiddles(N2, N1);
	std::array<complex<fft_simd4>, N1> p;
	for (int k2 = 0; k2 < N2 / 2; k2 += SIMD_SIZE) {
		for (int si = 0; si < SIMD_SIZE; si++) {
			for (int n1 = 0; n1 < N1; n1++) {
				p[n1].real()[si] = y[N1 * (k2 + si) + n1];
				p[n1].imag()[si] = y[N - N1 * (k2 + si) + n1];
			}
		}
		for (int n1 = 1; n1 < N1; n1++) {
			p[n1] *= W[k2 / SIMD_SIZE][n1];
		}
		sfft_complex<N1>((fft_simd4*) p.data());
		for (int k1 = 0; k1 < N1 / 2; k1++) {
			const int k = N2 * k1 + k2;
			for (int si = 0; si < SIMD_SIZE; si++) {
				x[k + si] = p[k1].real()[si];
				x[N - k - si] = p[k1].imag()[si];
			}
		}
		for (int k1 = N1 / 2; k1 < N1; k1++) {
			const int k = N2 * k1 + k2;
			for (int si = 0; si < SIMD_SIZE; si++) {
				x[N - k - si] = p[k1].real()[si];
				x[k + si] = -p[k1].imag()[si];
			}
		}
	}
	std::array<double, N1> q;
	for (int n1 = 0; n1 < N1; n1++) {
		q[n1] = y[n1];
	}
	sfft_real<N1>(q.data());
	x[0] = q[0];
	x[N / 2] = q[N1 / 2];
	for (int k1 = 1; k1 < N1 / 2; k1++) {
		x[N2 * k1] = q[k1];
		x[N - N2 * k1] = q[N1 - k1];
	}
	for (int n1 = 0; n1 < N1; n1++) {
		q[n1] = y[N / 2 + n1];
	}
	sfft_skew<N1>(q.data());
	for (int k1 = 0; k1 < N1 / 2; k1++) {
		x[(N2 * k1 + N2 / 2)] = q[k1];
		x[N - (N2 * k1 + N2 / 2)] = q[N1 - k1 - 1];
	}
}

