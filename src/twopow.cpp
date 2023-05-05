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
void fft3(T* X, int N) {
	const auto& W = twiddles(N);
	constexpr int M1 = 2;
	static std::vector<int> D;
	if (!D.size()) {
		D.resize(M1);
		std::iota(D.begin(), D.end(), 0);
		scramble(D.data(), M1);
	}
	T z[2 * M1];
	T* q = (T*) z;
	complex<T>* p = (complex<T>*) z;
	int M2 = 1;
	while (M2 < N) {
		const int M = M1 * M2;
		const int L = N / M;
		for (int l = 0; l < L; l++) {
			const int L0 = l * M;
			for (int n1 = 0; n1 < M1; n1++) {
				q[n1] = X[L0 + M2 * D[n1]];
			}
			sfft_real<M1>(q);
			X[L0 + 0] = q[0];
			for (int k1 = 1; k1 < M1 / 2; k1++) {
				X[L0 + 0 + M2 * k1] = q[k1];
				X[L0 + M - M2 * k1] = q[M1 - k1];
			}
			X[L0 + M / 2] = q[M1 / 2];
			if (M2 / 2) {
				for (int n1 = 0; n1 < M1; n1++) {
					q[n1] = X[L0 + M2 * D[n1] + M2 / 2];
				}
				sfft_skew<M1>(q);
				for (int k1 = 0; k1 < M1 / 2; k1++) {
					X[L0 + 0 + (M2 * k1 + M2 / 2)] = q[k1];
					X[L0 + M - (M2 * k1 + M2 / 2)] = q[M1 - k1 - 1];
				}
			}
			for (int n1 = 1; n1 < M1; n1++) {
				const auto m1 = D[n1];
				const int M2m1 = M2 * m1;
				const int M2m1pM2 = M2m1 + M2;
				for (int k2 = 1; k2 < M2 / 2; k2++) {
					auto& re = X[L0 + M2m1 + k2];
					auto& im = X[L0 + M2m1pM2 - k2];
					const auto re0 = re;
					const auto c = W[n1 * k2 * L].real();
					const auto s = W[n1 * k2 * L].imag();
					re = re * c - im * s;
					im = re0 * s + im * c;
				}
			}
			for (int k2 = 1; k2 < M2 / 2; k2++) {
				for (int n1 = 0; n1 < M1; n1++) {
					p[n1].real() = X[L0 + M2 * D[n1] + k2];
					p[n1].imag() = X[L0 + M2 * D[n1] - k2 + M2];
				}
				sfft_complex<M1>((T*) p);
				for (int k1 = 0; k1 < M1 / 2; k1++) {
					const int k = M2 * k1 + k2;
					X[L0 + k] = p[k1].real();
					X[L0 + M - k] = p[k1].imag();
				}
				for (int k1 = M1 / 2; k1 < M1; k1++) {
					const int k = M2 * k1 + k2;
					X[L0 + M - k] = p[k1].real();
					X[L0 + k] = -p[k1].imag();
				}
			}
		}
		M2 *= M1;
	}
}

template<class T>
void fft2(T* X, int N) {
	constexpr int N1 = 4;
	const int N2 = N / N1;

	T T0R, T0I, T1R, T1I, T2R, T2I, T3R, T3I;
	T U0R, U0I, U1R, U1I, U2R, U2I, U3R, U3I;
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
		fft2(&X[n1 * N2], N2);
	}

	const auto& W = twiddles(N);

	i0 = 0;
	i1 = i0 + N2;
	i2 = i1 + N2;
	i3 = i2 + N2;
	U0R = X[i0];
	U2R = X[i1];
	U1R = X[i2];
	U3R = X[i3];
	T0R = U0R + U2R;
	T2R = U0R - U2R;
	T1R = U1R + U3R;
	T3R = U1R - U3R;
	X[i0] = T0R + T1R;
	X[i1] = T2R;
	X[i2] = T0R - T1R;
	X[i3] = -T3R;
	if (N >= 8) {
		i0 = N2 / 2;
		i1 = i0 + N2;
		i2 = i1 + N2;
		i3 = i2 + N2;
		U0R = X[i0];
		U2R = X[i1];
		U1R = X[i2];
		U3R = X[i3];
		T1R = (U1R - U3R) * M_SQRT1_2;
		T3R = (U1R + U3R) * M_SQRT1_2;
		X[i0] = U0R + T1R;
		X[i1] = U0R - T1R;
		X[i2] = U2R - T3R;
		X[i3] = -U2R - T3R;
	}
	for (int k2 = 1; k2 < N2 / 2; k2++) {
		i0 = k2;
		i1 = i0 + N2;
		i2 = i1 + N2;
		i3 = i2 + N2;
		i4 = N2 - k2;
		i5 = i4 + N2;
		i6 = i5 + N2;
		i7 = i6 + N2;
		U0R = X[i0];
		U2R = X[i1];
		U1R = X[i2];
		U3R = X[i3];
		U0I = X[i4];
		U1I = X[i6];
		U2I = X[i5];
		U3I = X[i7];
		j1 = k2;
		j2 = 2 * k2;
		j3 = 3 * k2;
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
		T0R = U0R + U2R;
		T0I = U0I + U2I;
		T2R = U0R - U2R;
		T2I = U0I - U2I;
		T1R = U1R + U3R;
		T1I = U1I + U3I;
		T3R = U1R - U3R;
		T3I = U1I - U3I;
		X[i0] = T0R + T1R;
		X[i1] = T2R + T3I;
		X[i2] = -T0I + T1I;
		X[i3] = -T2I - T3R;
		X[i4] = T2R - T3I;
		X[i5] = T0R - T1R;
		X[i6] = T2I - T3R;
		X[i7] = T0I + T1I;
	}
}

void fft_2pow(double* x, int N) {
	scramble((double*) x, N);
	fft2((double*) x, N);
}

