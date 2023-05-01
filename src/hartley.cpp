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
	y.resize(N / 2);
	double t1 = 0.0;
	double t2 = 0.0;
	const auto& w = twiddles(2 * N);
	for (int n = 0; n < N; n += 2) {
		t1 += x[n] * w[n].real();
	}
	for (int n = 1; n < N; n += 2) {
		t2 += x[n] * w[n].real();
	}
	for (int n = 0; n < N; n++) {
		x[n] *= 4.0 * w[n].imag();
	}
	fft_real(x, N);
	y[0].real() = 2.0 * (t1 + t2);
	y[0].imag() = -0.5 * x[0];
	for (int n = 1; n < N / 2 - 1; n++) {
		y[n].real() = -x[N - n] + y[n - 1].real();
		y[n].imag() = -x[n] + y[n - 1].imag();
	}
	y[N / 2 - 1].real() = 2.0 * (t1 - t2);
	y[N / 2 - 1].imag() = 0.5 * x[N / 2];
	for (int n = 0; n < N / 2; n++) {
		x[n] = 0.5 * y[n].real();
		x[N - n - 1] = -0.5 * y[n].imag();
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
	int L = ilogb(N);
	if (N != 1 << L || L % 2 == 1) {
		printf("6 step not eligible\n");
		abort();
	}
	const auto& W = twiddles(N);
	const int M = 1 << (L / 2);
	const int M2 = M / 2;
	const int N2 = N / 2;
	for (int n1 = 0; n1 < M; n1++) {
		for (int n2 = n1 + 1; n2 < M; n2++) {
			std::swap(X[M * n2 + n1], X[M * n1 + n2]);
		}
	}

	for (int n1 = 0; n1 < M; n1++) {
		fft_real(X + M * n1, M);
	}

	for (int k2 = 1; k2 < M2; k2++) {
		for (int n1 = 1; n1 < M; n1++) {
			const auto w = W[n1 * k2];
			const auto r = M * n1 + k2;
			const auto i = M * n1 + (M - k2);
			auto xr = X[r];
			auto xi = X[i];
			X[r] = xr * w.real() - xi * w.imag();
			X[i] = xr * w.imag() + xi * w.real();
		}
	}

	for (int k2 = 0; k2 < M; k2++) {
		for (int n1 = k2 + 1; n1 < M; n1++) {
			std::swap(X[M * n1 + k2], X[M * k2 + n1]);
		}
	}

	for (int k2 = 0; k2 < M; k2++) {
		if (k2 != M / 2) {
			fft_real(X + M * k2, M);
		} else {
			fft_real_odd(X + M * k2, M);
		}
	}

	for (int k2 = 1; k2 < M2; k2++) {
		for (int k1 = 1; k1 < M2; k1++) {
			const int rr = M * k2 + k1;
			const int ri = M * k2 + M - k1;
			const int ir = N - M * k2 + k1;
			const int ii = N - M * k2 + M - k1;
			const auto xrr = X[rr];
			const auto xri = X[ri];
			const auto xir = X[ir];
			const auto xii = X[ii];
			X[rr] = xrr - xii;
			X[ri] = xrr + xii;
			X[ir] = xir + xri;
			X[ii] = xir - xri;

		}
	}
	for (int k2 = 0; k2 < M2; k2++) {
		for (int k1 = 0; k1 < M2; k1++) {
			std::swap(X[N2 + M * k2 + (M - k1) % M], X[N2 + M2 + M * k2 + k1]);
			X[N2 + M * k2 + (M - k1) % M] = -X[N2 + M * k2 + (M - k1) % M];
		}
	}
	for (int k2 = 0; k2 < M2; k2++) {
		for (int k1 = k2 + 1; k1 < M2; k1++) {
			const int i1 = M * k2 + k1;
			const int i2 = M * k1 + k2;
			std::swap(X[i1], X[i2]);
			std::swap(X[i1 + M2], X[i2 + M2]);
			std::swap(X[i1 + N2], X[i2 + N2]);
			std::swap(X[i1 + M2 + N2], X[i2 + M2 + N2]);
		}
	}
	for (int k2 = 1; k2 < M2; k2++) {
		for (int k1 = 0; k1 < M2; k1++) {
			int i = M * k1 + k2;
			int j = M * ((M2 - k1) % M2) + (M2 - k2);
			if (i < j) {
				std::swap(X[i + M2], X[j + M2]);
			}
		}
	}
	for (int k2 = 1; k2 < M2; k2++) {
		for (int k1 = 0; k1 < M2; k1++) {
			int i = M * k1 + k2;
			int j = M * ((M2 - k1) % M2) + (M2 - k2) % M2;
			if (i < j) {
				std::swap(X[i + N2], X[j + N2]);
			}
		}
	}
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

void fft_simd_8(__m256d* x, __m256d* y) {

	static constexpr __m256d TC1 = {1.0,cos(1.0 * 2.0*M_PI/8.0),cos(2.0 * 2.0*M_PI/8.0),cos(3.0 * 2.0*M_PI/8.0)};
	static constexpr __m256d TS1 = {0.0,sin(-1.0 * 2.0*M_PI/8.0),sin(-2.0 * 2.0*M_PI/8.0),sin(-3.0 * 2.0*M_PI/8.0)};
	static constexpr int scl = 8;
	static constexpr __m256i I = {0, 2, 4, 6};
	static constexpr __m256d T = {1.0, -1.0, 1.0, -1.0};

	auto& R0 = x[0];
	auto& R1 = x[1];
	auto& I0 = y[0];
	auto& I1 = y[1];

	__m256d A0, A1;
	__m256d B0, B1;

	A0 = _mm256_i64gather_pd((((double*)x) + 0), I, scl);
	A1 = _mm256_i64gather_pd((((double*)x) + 1), I, scl);
	B0 = _mm256_i64gather_pd((((double*)y) + 0), I, scl);
	B1 = _mm256_i64gather_pd((((double*)y) + 1), I, scl);
	R0 = _mm256_add_pd(A0, A1);
	R1 = _mm256_sub_pd(A0, A1);
	I0 = _mm256_add_pd(B0, B1);
	I1 = _mm256_sub_pd(B0, B1);
	std::swap(R1[1], I1[1]);
	std::swap(R1[3], I1[3]);
	I1 = _mm256_mul_pd(I1, T);
	A0 = _mm256_i64gather_pd((((double*)x) + 0), I, scl);
	A1 = _mm256_i64gather_pd((((double*)x) + 1), I, scl);
	B0 = _mm256_i64gather_pd((((double*)y) + 0), I, scl);
	B1 = _mm256_i64gather_pd((((double*)y) + 1), I, scl);

	R0 = _mm256_add_pd(A0, A1);
	R1 = _mm256_sub_pd(A0, A1);
	I0 = _mm256_add_pd(B0, B1);
	I1 = _mm256_sub_pd(B0, B1);
	A0 = _mm256_i64gather_pd((((double*)x) + 0), I, scl);
	A1 = _mm256_i64gather_pd((((double*)x) + 1), I, scl);
	B0 = _mm256_i64gather_pd((((double*)y) + 0), I, scl);
	B1 = _mm256_i64gather_pd((((double*)y) + 1), I, scl);
	R0 = A0;
	R1 = A1;
	I0 = B0;
	I1 = B1;
	A1 = _mm256_mul_pd(I1, TS1);
	A1 = _mm256_fmsub_pd(R1, TC1, A1);
	B1 = _mm256_mul_pd(I1, TC1);
	B1 = _mm256_fmadd_pd(R1, TS1, B1);
	R0 = _mm256_add_pd(A0, A1);
	R1 = _mm256_sub_pd(A0, A1);
	I0 = _mm256_add_pd(B0, B1);
	I1 = _mm256_sub_pd(B0, B1);
}

void fft_simd_16(__m256d* x, __m256d* y) {

	static constexpr __m256d TC1 = {1.0,cos(1.0 * 2.0*M_PI/16.0),cos(2.0 * 2.0*M_PI/16.0),cos(3.0 * 2.0*M_PI/16.0)};
	static constexpr __m256d TC2 = {1.0,cos(2.0 * 2.0*M_PI/16.0),cos(4.0 * 2.0*M_PI/16.0),cos(6.0 * 2.0*M_PI/16.0)};
	static constexpr __m256d TC3 = {1.0,cos(3.0 * 2.0*M_PI/16.0),cos(6.0 * 2.0*M_PI/16.0),cos(9.0 * 2.0*M_PI/16.0)};
	static constexpr __m256d TS1 = {0.0,-sin(1.0 * 2.0*M_PI/16.0),-sin(2.0 * 2.0*M_PI/16.0),-sin(3.0 * 2.0*M_PI/16.0)};
	static constexpr __m256d TS2 = {0.0,-sin(2.0 * 2.0*M_PI/16.0),-sin(4.0 * 2.0*M_PI/16.0),-sin(6.0 * 2.0*M_PI/16.0)};
	static constexpr __m256d TS3 = {0.0,-sin(3.0 * 2.0*M_PI/16.0),-sin(6.0 * 2.0*M_PI/16.0),-sin(9.0 * 2.0*M_PI/16.0)};
	static constexpr int scl = 8;
	static constexpr __m256i I = {0, 4, 8, 12};

	auto& R0 = x[0];
	auto& R1 = x[1];
	auto& R2 = x[2];
	auto& R3 = x[3];
	auto& I0 = y[0];
	auto& I1 = y[1];
	auto& I2 = y[2];
	auto& I3 = y[3];

	__m256d A0, A1, A2, A3;
	__m256d B0, B1, B2, B3;

	A0 = _mm256_i64gather_pd((((double*)x) + 0), I, scl);
	A1 = _mm256_i64gather_pd((((double*)x) + 1), I, scl);
	A2 = _mm256_i64gather_pd((((double*)x) + 2), I, scl);
	A3 = _mm256_i64gather_pd((((double*)x) + 3), I, scl);
	B0 = _mm256_i64gather_pd((((double*)y) + 0), I, scl);
	B1 = _mm256_i64gather_pd((((double*)y) + 1), I, scl);
	B2 = _mm256_i64gather_pd((((double*)y) + 2), I, scl);
	B3 = _mm256_i64gather_pd((((double*)y) + 3), I, scl);
	R0 = A0;
	R1 = A1;
	R2 = A2;
	R3 = A3;
	I0 = B0;
	I1 = B1;
	I2 = B2;
	I3 = B3;
	A0 = _mm256_add_pd(R0, R1);
	A1 = _mm256_sub_pd(R0, R1);
	A2 = _mm256_add_pd(R2, R3);
	A3 = _mm256_sub_pd(R2, R3);
	B0 = _mm256_add_pd(I0, I1);
	B1 = _mm256_sub_pd(I0, I1);
	B2 = _mm256_add_pd(I2, I3);
	B3 = _mm256_sub_pd(I2, I3);
	R0 = _mm256_add_pd(A0, A2);
	R2 = _mm256_sub_pd(A0, A2);
	R1 = _mm256_add_pd(A1, B3);
	R3 = _mm256_sub_pd(A1, B3);
	I0 = _mm256_add_pd(B0, B2);
	I2 = _mm256_sub_pd(B0, B2);
	I1 = _mm256_sub_pd(B1, A3);
	I3 = _mm256_add_pd(B1, A3);
	A0 = _mm256_i64gather_pd((((double*)x) + 0), I, scl);
	A1 = _mm256_i64gather_pd((((double*)x) + 1), I, scl);
	A2 = _mm256_i64gather_pd((((double*)x) + 2), I, scl);
	A3 = _mm256_i64gather_pd((((double*)x) + 3), I, scl);
	B0 = _mm256_i64gather_pd((((double*)y) + 0), I, scl);
	B1 = _mm256_i64gather_pd((((double*)y) + 1), I, scl);
	B2 = _mm256_i64gather_pd((((double*)y) + 2), I, scl);
	B3 = _mm256_i64gather_pd((((double*)y) + 3), I, scl);
	R0 = A0;
	R1 = _mm256_mul_pd(B1, TS2);
	R2 = _mm256_mul_pd(B2, TS1);
	R3 = _mm256_mul_pd(B3, TS3);
	R1 = _mm256_fmsub_pd(A1, TC2, R1);
	R2 = _mm256_fmsub_pd(A2, TC1, R2);
	R3 = _mm256_fmsub_pd(A3, TC3, R3);
	I0 = B0;
	I1 = _mm256_mul_pd(A1, TS2);
	I2 = _mm256_mul_pd(A2, TS1);
	I3 = _mm256_mul_pd(A3, TS3);
	I1 = _mm256_fmadd_pd(B1, TC2, I1);
	I2 = _mm256_fmadd_pd(B2, TC1, I2);
	I3 = _mm256_fmadd_pd(B3, TC3, I3);
	A0 = _mm256_add_pd(R0, R1);
	A1 = _mm256_sub_pd(R0, R1);
	A2 = _mm256_add_pd(R2, R3);
	A3 = _mm256_sub_pd(R2, R3);
	B0 = _mm256_add_pd(I0, I1);
	B1 = _mm256_sub_pd(I0, I1);
	B2 = _mm256_add_pd(I2, I3);
	B3 = _mm256_sub_pd(I2, I3);
	R0 = _mm256_add_pd(A0, A2);
	R2 = _mm256_sub_pd(A0, A2);
	R1 = _mm256_add_pd(A1, B3);
	R3 = _mm256_sub_pd(A1, B3);
	I0 = _mm256_add_pd(B0, B2);
	I2 = _mm256_sub_pd(B0, B2);
	I1 = _mm256_sub_pd(B1, A3);
	I3 = _mm256_add_pd(B1, A3);

}

constexpr double C(int n, int N) {
	return cos(2.0 * M_PI * n / N);
}

constexpr double S(int n, int N) {
	return -sin(2.0 * M_PI * n / N);
}

void fft_simd_32(__m256d* x, __m256d* y) {

	static constexpr __m256d TC1 = {1.0,C(1,32), C(2,32), C(3,32)};
	static constexpr __m256d TS1 = {0.0,S(1,32), S(2,32), S(3,32)};
	static constexpr __m256d TC2 = {C(4,32), C(5,32), C(6,32), C(7,32)};
	static constexpr __m256d TS2 = {S(4,32), S(5,32), S(6,32), S(7,32)};
	static constexpr __m256d TC3 = {C(8,32), C(9,32), C(10,32), C(11,32)};
	static constexpr __m256d TS3 = {S(8,32), S(9,32), S(10,32), S(11,32)};
	static constexpr __m256d TC4 = {C(12,32), C(13,32), C(14,32), C(15,32)};
	static constexpr __m256d TS4 = {S(12,32), S(13,32), S(14,32), S(15,32)};

	auto& R0 = x[0];
	auto& R1 = x[1];
	auto& R2 = x[2];
	auto& R3 = x[3];
	auto& R4 = x[4];
	auto& R5 = x[5];
	auto& R6 = x[6];
	auto& R7 = x[7];
	auto& I0 = y[0];
	auto& I1 = y[1];
	auto& I2 = y[2];
	auto& I3 = y[3];
	auto& I4 = y[4];
	auto& I5 = y[5];
	auto& I6 = y[6];
	auto& I7 = y[7];

	__m256d A0, A1, A2, A3, A4, A5, A6, A7;
	__m256d B0, B1, B2, B3, B4, B5, B6, B7;

	fft_simd_16(x, y);
	fft_simd_16(x + 4,y + 4);

	A0 = R0;
	A1 = R1;
	A2 = R2;
	A3 = R3;
	A4 = _mm256_mul_pd(I4, TS1);
	A4 = _mm256_fmsub_pd(R4, TC1, A4);
	A5 = _mm256_mul_pd(I5, TS2);
	A5 = _mm256_fmsub_pd(R5, TC2, A5);
	A6 = _mm256_mul_pd(I6, TS3);
	A6 = _mm256_fmsub_pd(R6, TC3, A6);
	A7 = _mm256_mul_pd(I7, TS4);
	A7 = _mm256_fmsub_pd(R7, TC4, A7);

	B0 = I0;
	B1 = I1;
	B2 = I2;
	B3 = I3;
	B4 = _mm256_mul_pd(R4, TS1);
	B4 = _mm256_fmadd_pd(I4, TC1, B4);
	B5 = _mm256_mul_pd(R5, TS2);
	B5 = _mm256_fmadd_pd(I5, TC2, B5);
	B6 = _mm256_mul_pd(R6, TS3);
	B6 = _mm256_fmadd_pd(I6, TC3, B6);
	B7 = _mm256_mul_pd(R7, TS4);
	B7 = _mm256_fmadd_pd(I7, TC4, B7);

	R0 = _mm256_add_pd(A0, A4);
	R4 = _mm256_sub_pd(A0, A4);
	R1 = _mm256_add_pd(A1, A5);
	R5 = _mm256_sub_pd(A1, A5);
	R2 = _mm256_add_pd(A2, A6);
	R6 = _mm256_sub_pd(A2, A6);
	R3 = _mm256_add_pd(A3, A7);
	R7 = _mm256_sub_pd(A3, A7);

	I0 = _mm256_add_pd(B0, B4);
	I4 = _mm256_sub_pd(B0, B4);
	I1 = _mm256_add_pd(B1, B5);
	I5 = _mm256_sub_pd(B1, B5);
	I2 = _mm256_add_pd(B2, B6);
	I6 = _mm256_sub_pd(B2, B6);
	I3 = _mm256_add_pd(B3, B7);
	I7 = _mm256_sub_pd(B3, B7);

}

#include <memory>

const std::vector<std::vector<complex<__m256d >>>& simd_twiddles(int N1, int N2) {
	using entry_type = std::shared_ptr<std::vector<std::vector<complex<__m256d>>>>;
	static std::unordered_map<int, std::unordered_map<int, entry_type>> cache;
	auto iter = cache[N1].find(N2);
	if (iter != cache[N1].end()) {
		return *(iter->second);
	} else {
		const int N = N1 * N2;
		const int N1v = round_down(N1, SIMD_SIZE);
		std::vector<std::vector<complex<__m256d>>>W(N1v/SIMD_SIZE, std::vector<complex<__m256d>>(N2));
		for (int n = 0; n < N1v; n += SIMD_SIZE) {
			for (int k = 0; k < N2; k++) {
				for (int i = 0; i < SIMD_SIZE; i++) {
					W[n / SIMD_SIZE][k].real()[i] = cos(-2.0 * M_PI * (n+i) * k / N);
					W[n / SIMD_SIZE][k].imag()[i] = sin(-2.0 * M_PI * (n+i) * k / N);
				}
			}
		}
		cache[N1][N2] = std::make_shared<std::vector<std::vector<complex<__m256d>>>>(std::move(W));
		return *(cache[N1][N2]);
	}
}

void fft_2pow(double* x, double* y, int N) {
	constexpr int N1 = 8;
	if (N == 8) {
		fft_simd_8((__m256d *) x, (__m256d *) y);
		return;
	} else if (N == 16) {
		fft_simd_16((__m256d *) x, (__m256d *) y);
		return;
	} else if (N == 32) {
		fft_simd_32((__m256d *) x, (__m256d *) y);
		return;
	}
	constexpr int bit_reverse[N1] = { 0, 4, 2, 6, 1, 5, 3, 7 };
	const int N2 = N / N1;
	const auto& W = simd_twiddles(N1, N2);
	for (int n1 = 0; n1 < N1; n1++) {
		fft_2pow(x + N2 * n1, y + N2 * n1, N2);
	}
	double re[N1];
	double im[N1];
	for (int k2 = 0; k2 < N2; k2++) {
		for (int n1 = 0; n1 < N1; n1++) {
			const int i = N2 * n1 + k2;
			re[n1] = x[i];
			im[n1] = y[i];
		}
		for (int n1 = 0; n1 < N1; n1++) {
			int j = bit_reverse[n1];
			if( n1 < j ) {
				std::swap(re[n1], re[j]);
				std::swap(im[n1], im[j]);
			}
		}
		for (int n1 = 0, n0 = 0; n1 < N1; n1 += SIMD_SIZE, n0++) {
			__m256d xre = *((__m256d *) &re[n1]);
			__m256d xim = *((__m256d *) &im[n1]);
			const __m256d& wre = W[n0][k2].real();
			const __m256d& wim = W[n0][k2].imag();
			*((__m256d *) &re[n1]) = _mm256_fmsub_pd(xre, wre, _mm256_mul_pd(xim, wim));
			*((__m256d *) &im[n1]) = _mm256_fmadd_pd(xim, wre, _mm256_mul_pd(xre, wim));
		}
		for (int n1 = 0; n1 < N1; n1++) {
			int j = bit_reverse[n1];
			if( n1 < j ) {
				std::swap(re[n1], re[j]);
				std::swap(im[n1], im[j]);
			}
		}
		switch (N1) {
		case 8:
			fft_simd_8((__m256d *) re, (__m256d *) im);
			break;
		case 16:
			fft_simd_16((__m256d *) re, (__m256d *) im);
			break;
		case 32:
			fft_simd_32((__m256d *) re, (__m256d *) im);
			break;
		}
		for (int n1 = 0; n1 < N1; n1++) {
			const int i = N2 * n1 + k2;
			x[i] = re[n1];
			y[i] = im[n1];
		}
	}
}

template<class T>
void fft_split_real_2pow(T* X, int N) {

	if (N <= 64) {
		sfft_real_scrambled(X, N);
		return;
	}

	const int No2 = N >> 1;
	const int No4 = N >> 2;
	const int No8 = N >> 3;
	const auto& W = twiddles(N);

	fft_split_real_2pow(X, No2);
	fft_split_real_2pow(X + No2, No4);
	fft_split_real_2pow(X + No2 + No4, No4);

	int I0, J0, J1, J2, J3;
	int I1 = No4;
	int I2 = I1 + No4;
	int I3 = I2 + No4;
	auto qe0 = X[0];
	auto qe1 = X[I1];
	auto qo0 = X[I2];
	auto qo1 = X[I3];
	auto t0 = qo1 + qo0;
	auto t1 = qo1 - qo0;
	X[0] = qe0 + t0;
	X[I1] = qe1;
	X[I2] = qe0 - t0;
	X[I3] = t1;
	if (No8) {
		constexpr auto c0 = (-1.0 / M_SQRT2);
		I0 = No8;
		I1 = I0 + No4;
		I2 = I1 + No4;
		I3 = I2 + No4;
		qe0 = X[I0];
		qe1 = X[I1];
		qo0 = X[I2];
		qo1 = X[I3];
		t0 = (qo1 - qo0) * c0;
		t1 = (qo0 + qo1) * c0;
		X[I0] = qe0 + t0;
		X[I3] = qe1 + t1;
		X[I1] = qe0 - t0;
		X[I2] = t1 - qe1;
	}
	for (int k2 = 1; k2 < No8; k2++) {
		const auto& W1 = W[k2];
		const auto& W3 = W[3 * k2];
		const auto& C1 = W1.real();
		const auto& C3 = W3.real();
		const auto& S1 = W1.imag();
		const auto& S3 = W3.imag();
		I0 = k2;
		J0 = No4 - k2;
		I1 = I0 + No4;
		I2 = I1 + No4;
		I3 = I2 + No4;
		J1 = J0 + No4;
		J2 = J1 + No4;
		J3 = J2 + No4;
		auto qe0r = X[I0];
		auto qe0i = X[J1];
		auto qe1r = X[J0];
		auto qe1i = X[I1];
		auto qo0r = X[I2];
		auto qo0i = X[J2];
		auto qo1r = X[I3];
		auto qo1i = X[J3];
		t0 = qo0r;
		t1 = qo1r;
		qo0r = qo0r * C1 - qo0i * S1;
		qo1r = qo1r * C3 - qo1i * S3;
		qo0i = t0 * S1 + qo0i * C1;
		qo1i = t1 * S3 + qo1i * C3;
		const auto zsr = qo0r + qo1r;
		const auto zsi = qo0i + qo1i;
		const auto zdr = qo1r - qo0r;
		const auto zdi = qo1i - qo0i;
		X[I0] = qe0r + zsr;
		X[J3] = qe0i + zsi;
		X[I1] = qe1r - zdi;
		X[J2] = zdr - qe1i;
		X[J1] = qe0r - zsr;
		X[I2] = zsi - qe0i;
		X[J0] = qe1r + zdi;
		X[I3] = qe1i + zdr;
	}
}

