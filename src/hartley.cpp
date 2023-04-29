#include "fftu.hpp"
#include "util.hpp"
#include <numeric>

#include <cmath>
#include <numeric>
/*
 template<class T>
 void fht_radix2(T* X, int N) {
 if (N == 2) {
 const auto T1 = X[0];
 const auto T2 = X[1];
 X[0] = T1 + T2;
 X[1] = T1 - T2;
 return;
 } else if (N <= 1) {
 return;
 }
 const auto& C = cos_twiddles(N);
 const auto& S = cos_twiddles(N);
 const int N2 = N >> 1;
 const int N4 = N >> 2;
 const int N8 = N >> 3;
 const int L12 = N4;
 const int L13 = L12 + N4;
 const int L14 = L13 + N4;
 auto T1 = X[L12] - X[L14];
 auto T2 = X[0] + X[L13];
 auto T3 = X[L12] + X[L14];
 auto T4 = X[0] - X[L13];
 X[L14] = T4 - T1;
 X[L13] = T4 + T1;
 X[L12] = T3;
 X[0] = T2;
 if (N4 != 1) {
 const int L21 = N8;
 const int L22 = L21 + N4;
 const int L23 = L22 + N4;
 const int L24 = L23 + N4;
 T1 = X[L22];
 T2 = X[L23];
 T3 = X[L24];
 X[L24] = (T1 - T3) * M_SQRT2;
 X[L23] = (X[L21] - T2) * M_SQRT2;
 X[L22] = (T1 + T3);
 X[L21] = (X[L21] + T2);
 }
 for (int k = 1; k < N8; k++) {
 const int L11 = k;
 const int L12 = L11 + N4;
 const int L13 = L12 + N4;
 const int L14 = L13 + N4;
 const int L21 = N4 - k;
 const int L22 = L21 + N4;
 const int L23 = L22 + N4;
 const int L24 = L23 + N4;
 const auto& C1 = C[k];
 const auto& S1 = S[k];
 const auto& C3 = C[3 * k];
 const auto& S3 = S[3 * k];
 T T1, T2, T3, T4, T5;
 T5 = X[L21] - X[L23];
 T2 = X[L11] - X[L13];
 T1 = T2 + T5;
 T2 -= T5;
 T5 = X[L22] - X[L24];
 T4 = X[L14] - X[L12];
 T3 = T4 + T5;
 T4 -= T5;
 X[L11] += X[L13];
 X[L12] += X[L14];
 X[L21] += X[L23];
 X[L22] += X[L24];
 X[L13] = T1 * C1 + T3 * S1;
 X[L14] = T2 * C3 - T4 * S3;
 X[L23] = T1 * S1 - T3 * C1;
 X[L24] = T2 * S3 + T4 * C3;
 }
 fht_radix2<T>(X, N2);
 fht_radix2<T>(X + N2, N4);
 fht_radix2<T>(X + N2 + N4, N4);
 }
 */

template<class T>
void fht_radix2(T* x, int N, bool scramble = false) {
	if (scramble) {
		int j = 0;
		for (int i = 0; i < N - 1; i++) {
			if (i > j) {
				const auto tmp = x[i];
				x[i] = x[j];
				x[j] = tmp;
			}
			int k = N / 2;
			while (k <= j) {
				j -= k;
				k >>= 1;
			}
			j += k;
		}
	}
	if (N <= 16) {
		sfht(x, N);
		return;
	}
	const int No4 = N >> 2;
	const int No8 = N >> 3;
	fht_radix2<T>(x, No4);
	fht_radix2<T>(x + No4, No4);
	fht_radix2<T>(x + 2 * No4, No4);
	fht_radix2<T>(x + 3 * No4, No4);
	const auto& C = cos_twiddles(N);
	const auto& S = sin_twiddles(N);
	const int L11 = 0;
	const int L12 = L11 + No4;
	const int L13 = L12 + No4;
	const int L14 = L13 + No4;
	const int L21 = No8;
	const int L22 = L21 + No4;
	const int L23 = L22 + No4;
	const int L24 = L23 + No4;
	auto T1 = x[L11] + x[L13];
	auto T2 = x[L11] - x[L13];
	auto T3 = x[L12] + x[L14];
	auto T4 = x[L12] - x[L14];
	x[L11] = T1 + T3;
	x[L12] = T1 - T3;
	x[L13] = T2 + T4;
	x[L14] = T2 - T4;
	T1 = x[L21];
	T2 = x[L23] * M_SQRT2;
	T3 = x[L22];
	T4 = x[L24] * M_SQRT2;
	x[L21] = T1 + T2 + T3;
	x[L22] = T1 - T3 + T4;
	x[L23] = T1 - T2 + T3;
	x[L24] = T1 - T3 - T4;
	for (int k = 1; k < No8; k++) {
		const auto C1 = C[k];
		const auto S1 = S[k];
		const auto C2 = C[2 * k];
		const auto S2 = S[2 * k];
		const auto C3 = C[3 * k];
		const auto S3 = S[3 * k];
		const int L11 = k;
		const int L12 = L11 + No4;
		const int L13 = L12 + No4;
		const int L14 = L13 + No4;
		const int L21 = No4 - k;
		const int L22 = L21 + No4;
		const int L23 = L22 + No4;
		const int L24 = L23 + No4;
		const auto T12 = x[L13] * C1 + x[L23] * S1;
		const auto T13 = x[L12] * C2 + x[L22] * S2;
		const auto T14 = x[L14] * C3 + x[L24] * S3;
		const auto T22 = x[L13] * S1 - x[L23] * C1;
		const auto T23 = x[L12] * S2 - x[L22] * C2;
		const auto T24 = x[L14] * S3 - x[L24] * C3;
		auto T1 = x[L21] + T23;
		auto T2 = x[L21] - T23;
		auto T3 = T22 + T24;
		auto T4 = T12 - T14;
		x[L24] = T2 - T3;
		x[L23] = T1 - T4;
		x[L22] = T2 + T3;
		x[L21] = T1 + T4;
		T1 = x[L11] + T13;
		T2 = x[L11] - T13;
		T3 = T24 - T22;
		T4 = T12 + T14;
		x[L14] = T2 - T3;
		x[L13] = T1 - T4;
		x[L12] = T2 + T3;
		x[L11] = T1 + T4;
	}
}

#include <cstring>

void fft_real_odd(double* x, int N) {
	static std::vector<complex<double>> y;
	y.resize(N / 4);
	double t1 = 0.0;
	double t2 = 0.0;
	const auto& w = twiddles(N);
	for (int n = 0; n < N / 2; n += 2) {
		t1 += x[n] * w[n].real();
	}
	for (int n = 1; n < N / 2; n += 2) {
		t2 += x[n] * w[n].real();
	}
	for (int n = 0; n < N / 2; n++) {
		x[n] *= 4.0 * w[n].imag();
	}
	fft_real(x, N / 2);
	y[0].real() = 2.0 * (t1 + t2);
	y[0].imag() = -0.5 * x[0];
	for (int n = 1; n < N / 4 - 1; n++) {
		y[n].real() = -x[N / 2 - n] + y[n - 1].real();
		y[n].imag() = -x[n] + y[n - 1].imag();
	}
	y[N / 4 - 1].real() = 2.0 * (t1 - t2);
	y[N / 4 - 1].imag() = 0.5 * x[N / 4];
	for (int n = 0; n < N / 4; n++) {
		x[n] = 0.5 * y[n].real();
		x[N / 2 - n - 1] = -0.5 * y[n].imag();
	}
}
void fht_odd(double* x, int N) {
	fft_real_odd(x, N);
	for (int n = 0; n < N / 2 - n - 1; n++) {
		const auto r = x[n];
		const auto i = x[N / 2 - n - 1];
		x[n] = r - i;
		x[N / 2 - n - 1] = r + i;
	}
}

void fft_radix6step(double* X, int N) {
	int M = ilogb(N);
	if (N != 1 << M || M % 2 == 1) {
		printf("6 step not eligible\n");
		abort();
	}
	const auto& W = twiddles(N);
	const int L = 1 << (M / 2);

	for (int n1 = 0; n1 < L; n1++) {
		for (int n2 = n1 + 1; n2 < L; n2++) {
			std::swap(X[L * n1 + n2], X[L * n2 + n1]);
		}
	}
	for (int n = 0; n < L; n++) {
		fht_radix2<double>(X + L * n, L, true);
	}
	for (int n = 0; n < L; n++) {
		for (int k = n + 1; k < L; k++) {
			std::swap(X[L * n + k], X[L * k + n]);
		}
	}
	for (int n = 0; n < L; n++) {
		for (int k = 1; k < L - k; k++) {
			const auto w = W[n * k];
			complex<double> z;
			z.real() = X[n + k * L];
			z.imag() = X[n + (L - k) * L];
			z *= w;
			X[n + k * L] = z.real();
			X[n + (L - k) * L] = z.imag();
		}
	}
	for (int n = 0; n < L; n++) {
		if (n != L / 2) {
			fht_radix2<double>(X + n * L, L, true);
		} else {
			fht_odd(X + N / 2, 2 * L);
		}
	}
	for (int k = 1; k < L - k; k++) {
		const auto R = 0.5 * (X[k] + X[L - k]);
		const auto I = -0.5 * (X[k] - X[L - k]);
		X[k] = R;
		X[L - k] = I;
	}
	for (int k = 0; k < L - k - 1; k++) {
		const auto e = X[N / 2 + k];
		const auto o = X[N / 2 + L - k - 1];
		X[N / 2 + k] = 0.5 * (e + o);
		X[N / 2 + L - k - 1] = -0.5 * (e - o);
	}
	for (int k1 = 1; k1 < L - k1; k1++) {
		const auto e = X[k1 * L];
		const auto o = X[(L - k1) * L];
		X[k1 * L] = 0.5 * (e + o);
		X[(L - k1) * L] = 0.5 * (o - e);
	}
	{
		int k2 = L / 2;
		for (int k1 = 1; k1 < L - k1; k1++) {
			const auto e = X[k2 + k1 * L];
			const auto o = X[k2 + (L - k1) * L];
			X[k2 + k1 * L] = 0.5 * (e + o);
			X[k2 + (L - k1) * L] = 0.5 * (o - e);
		}
	}
	for (int k2 = 1; k2 < L / 2; k2++) {
		for (int k1 = 1; k1 < L / 2; k1++) {
			const auto ER = 0.5 * X[L * k2 + k1];
			const auto EI = 0.5 * X[L * k2 + (L - k1) % L];
			const auto OR = -0.5 * X[L * (L - k2) + (L - k1) % L];
			const auto OI = 0.5 * X[L * (L - k2) + k1];
			X[L * k2 + k1] = EI + OI;
			X[N - L * k2 + k1] = -ER - OR;

			X[L * k2 + L - k1] = ER - OR;

			X[N - L * k2 + L - k1] = -EI + OI;
		}
	}

	for (int k1 = 0; k1 < L; k1++) {
		for (int k2 = k1 + 1; k2 < L; k2++) {
			std::swap(X[L * k1 + k2], X[L * k2 + k1]);
		}
	}
	for (int k1 = 0; k1 < L / 2; k1++) {
		for (int k2 = L / 2 + 1; k2 < L; k2++) {
			std::swap(X[L * k1 + k2], X[N - L * k1 - k2]);
		}
	}
	for (int k1 = L / 2; k1 < L; k1++) {
		for (int k2 = 1; k2 < L - k2 - 1; k2++) {
			X[L * k1 + L - k2] = -X[L * k1 + L - k2];
			std::swap(X[L * k1 + k2], X[L * k1 + L - k2]);
		}
	}
	for (int k = 1; k < N - k; k++) {
		const auto r = X[k];
		const auto i = X[N - k];
		X[k] = r - i;
		X[N - k] = r + i;
	}
	fht2fft(X, N);
}

void fht_radix2(double* X, int N) {
	constexpr int N1 = 4;
	using simd_t = fft_simd4;
	const int N2 = N / N1;
	int j = 0;
	simd_t* Z = (simd_t*) X;
	for (int i = 0; i < N2 - 1; i++) {
		if (i > j) {
			const auto tmp = Z[i];
			Z[i] = Z[j];
			Z[j] = tmp;
		}
		int k = N2 / 2;
		while (k <= j) {
			j -= k;
			k >>= 1;
		}
		j += k;
	}
	if (N <= 32) {
		sfht(X, N);
		return;
	}
	fht_radix2<simd_t>(Z, N2);
	static std::vector<double> Y;
	Y.resize(N);
	for (int k2 = 0; k2 < N2; k2++) {
		for (int n1 = 0; n1 < N1; n1++) {
			Y[N2 * n1 + k2] = X[N1 * k2 + n1];
		}
	}
	const auto& C = cos_twiddles_array(N1, N2);
	const auto& S = sin_twiddles_array(N1, N2);
	const int L11 = 0;
	const int L12 = L11 + N2;
	const int L13 = L12 + N2;
	const int L14 = L13 + N2;
	const int L21 = N2 / 2;
	const int L22 = L21 + N2;
	const int L23 = L22 + N2;
	const int L24 = L23 + N2;
	auto T1 = Y[L11] + Y[L12];
	auto T2 = Y[L11] - Y[L12];
	auto T3 = Y[L13] + Y[L14];
	auto T4 = Y[L13] - Y[L14];
	X[L11] = T1 + T3;
	X[L12] = T1 - T3;
	X[L13] = T2 + T4;
	X[L14] = T2 - T4;
	T1 = Y[L21];
	T2 = Y[L22] * M_SQRT2;
	T3 = Y[L23];
	T4 = Y[L24] * M_SQRT2;
	X[L21] = T1 + T2 + T3;
	X[L22] = T1 - T3 + T4;
	X[L23] = T1 - T2 + T3;
	X[L24] = T1 - T3 - T4;
	for (int k = 1; k < N1; k++) {
		const auto C1 = C[1][k];
		const auto S1 = S[1][k];
		const auto C2 = C[2][k];
		const auto S2 = S[2][k];
		const auto C3 = C[3][k];
		const auto S3 = S[3][k];
		const int L11 = k;
		const int L12 = L11 + N2;
		const int L13 = L12 + N2;
		const int L14 = L13 + N2;
		const int L21 = N2 - k;
		const int L22 = L21 + N2;
		const int L23 = L22 + N2;
		const int L24 = L23 + N2;
		const auto T12 = Y[L12] * C1 + Y[L22] * S1;
		const auto T13 = Y[L13] * C2 + Y[L23] * S2;
		const auto T14 = Y[L14] * C3 + Y[L24] * S3;
		const auto T22 = Y[L12] * S1 - Y[L22] * C1;
		const auto T23 = Y[L13] * S2 - Y[L23] * C2;
		const auto T24 = Y[L14] * S3 - Y[L24] * C3;
		auto T1 = Y[L21] + T23;
		auto T2 = Y[L21] - T23;
		auto T3 = T22 + T24;
		auto T4 = T12 - T14;
		X[L24] = T2 - T3;
		X[L23] = T1 - T4;
		X[L22] = T2 + T3;
		X[L21] = T1 + T4;
		T1 = Y[L11] + T13;
		T2 = Y[L11] - T13;
		T3 = T24 - T22;
		T4 = T12 + T14;
		X[L14] = T2 - T3;
		X[L13] = T1 - T4;
		X[L12] = T2 + T3;
		X[L11] = T1 + T4;
	}
	for (int k = N1; k < N2 / 2; k += N1) {
		const simd_t C1 = *((simd_t*) &C[1][k]);
		const simd_t S1 = *((simd_t*) &S[1][k]);
		const simd_t C2 = *((simd_t*) &C[2][k]);
		const simd_t S2 = *((simd_t*) &S[2][k]);
		const simd_t C3 = *((simd_t*) &C[3][k]);
		const simd_t S3 = *((simd_t*) &S[3][k]);
		const int L11 = k;
		const int L12 = L11 + N2;
		const int L13 = L12 + N2;
		const int L14 = L13 + N2;
		const int L21 = N2 - k;
		const int L22 = L21 + N2;
		const int L23 = L22 + N2;
		const int L24 = L23 + N2;
		const simd_t Y11 = *((simd_t*) &Y[L11]);
		const simd_t Y12 = *((simd_t*) &Y[L12]);
		const simd_t Y13 = *((simd_t*) &Y[L13]);
		const simd_t Y14 = *((simd_t*) &Y[L14]);
		simd_t Y21, Y22, Y23, Y24;
		for (int i = 0; i < N1; i++) {
			Y21[i] = Y[L21 - i];
			Y22[i] = Y[L22 - i];
			Y23[i] = Y[L23 - i];
			Y24[i] = Y[L24 - i];
		}
		const simd_t T12 = Y12 * C1 + Y22 * S1;
		const simd_t T13 = Y13 * C2 + Y23 * S2;
		const simd_t T14 = Y14 * C3 + Y24 * S3;
		const simd_t T22 = Y12 * S1 - Y22 * C1;
		const simd_t T23 = Y13 * S2 - Y23 * C2;
		const simd_t T24 = Y14 * S3 - Y24 * C3;
		simd_t T1 = Y21 + T23;
		simd_t T2 = Y21 - T23;
		simd_t T3 = T22 + T24;
		simd_t T4 = T12 - T14;
		const simd_t X24 = T2 - T3;
		const simd_t X23 = T1 - T4;
		const simd_t X22 = T2 + T3;
		const simd_t X21 = T1 + T4;
		T1 = Y11 + T13;
		T2 = Y11 - T13;
		T3 = T24 - T22;
		T4 = T12 + T14;
		const simd_t X14 = T2 - T3;
		const simd_t X13 = T1 - T4;
		const simd_t X12 = T2 + T3;
		const simd_t X11 = T1 + T4;
		for (int i = 0; i < N1; i++) {
			X[L21 - i] = X21[i];
			X[L22 - i] = X22[i];
			X[L23 - i] = X23[i];
			X[L24 - i] = X24[i];
		}
		*((simd_t*) &X[L11]) = X11;
		*((simd_t*) &X[L12]) = X12;
		*((simd_t*) &X[L13]) = X13;
		*((simd_t*) &X[L14]) = X14;
	}
}

