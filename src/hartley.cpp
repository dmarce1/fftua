#include "fftu.hpp"
#include "util.hpp"
#include <numeric>

#include <cmath>
#include <numeric>

void fht_radix2_indices(int* I, int N) {
	const int No4 = N >> 2;
	if (N <= SFFT_NMAX) {
		return;
	}
	std::vector<int> J(I, I + N);
	for (int n = 0; n < N / 4; n++) {
		I[n + 0 * N / 4] = J[4 * n + 0];
		I[n + 1 * N / 4] = J[4 * n + 1];
		I[n + 2 * N / 4] = J[4 * n + 2];
		I[n + 3 * N / 4] = J[4 * n + 3];
	}
	fht_radix2_indices(I, No4);
	fht_radix2_indices(I + No4, No4);
	fht_radix2_indices(I + 2 * No4, No4);
	fht_radix2_indices(I + 3 * No4, No4);
}

/*template<class T>
 void fht_radix2(T* x, int N, bool scramble = false) {
 const int No2 = N >> 1;
 if (scramble) {
 int j = 0;
 for (int i = 0; i < N - 1; i++) {
 if (i > j) {
 std::swap(x[i], x[j]);
 }
 int k = No2;
 while (k <= j) {
 j -= k;
 k >>= 1;
 }
 j += k;
 }
 }
 if (N == 64) {
 fht_64(x);
 return;
 } else if (N == 32) {
 fht_32(x);
 return;
 } else if (N < 32) {
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
 T2 = x[L23] * sqrt(2.0);
 T3 = x[L22];
 T4 = x[L24] * sqrt(2.0);
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
void fft_radix2(double* X, int N) {
	if (N == 2) {
		const auto z0 = X[0];
		const auto z1 = X[1];
		X[0] = z0 + z1;
		X[1] = z0 - z1;
		return;
	} else if (N <= 1) {
		return;
	}
	constexpr int N1 = 4;
	const int No2 = N / 2;
	const int No4 = N / 4;
	constexpr int N1o2 = N1 / 2;
	constexpr int N1o4 = N1 / 4;
	const int N2 = N / N1;
	const int N2o2 = N2 / 2;
	const auto& W = twiddles(N);
	fft_real(X, No2);
	for (int n1 = 0; n1 < N1o2; n1++) {
		int o = No2 + N2 * n1;
		fft_real(X + o, N2);
	}
	std::array<T, N1o2> qo;
	std::array<T, N1o2> qe;
	for (int n1 = 0; n1 < N1o2; n1++) {
		qe[n1] = X[N2 * n1];
	}
	for (int n1 = 0; n1 < N1o2; n1++) {
		qo[n1] = X[No2 + N2 * n1];
	}
	sfft_real_odd<N1>(qo.data());
	X[0] = qe[0] + qo[0];
	for (int k1 = 1; k1 < N1o4; k1++) {
		X[N2 * k1] = qe[k1] + qo[k1];
		X[N - N2 * k1] = qe[N1o2 - k1] + qo[N1o2 - k1];
	}
	X[No4] = qe[N1o4];
	X[N - No4] = qo[N1o4];
	for (int k1 = N1o4 + 1; k1 < N1o2; k1++) {
		X[N2 * k1] = qe[N1o2 - k1] - qo[N1o2 - k1];
		X[N - N2 * k1] = -qe[k1] + qo[k1];
	}
	X[No2] = qe[0] - qo[0];
	std::array<complex<T>, N1o2> po;
	std::array<complex<T>, N1o2> pe;
	for (int k2 = 1; k2 < (N2 + 1) / 2; k2++) {
		for (int n1 = 0; n1 < N1o2; n1++) {
			pe[n1].real() = X[N2 * n1 + k2];
			pe[n1].imag() = X[No2 - N2 * n1 - k2];
		}
		for (int n1 = 0; n1 < N1o2; n1++) {
			po[n1].real() = X[No2 + N2 * n1 + k2];
			po[n1].imag() = X[No2 + N2 * n1 - k2 + N2];
			po[n1] *= W[(2 * n1 + 1) * k2];
		}
		for (int n1 = (N1o2 + 1) / 2; n1 < N1o2; n1++) {
			std::swap(pe[n1].real(), pe[n1].imag());
			pe[n1] = pe[n1].conj();
		}
		sfft_complex_odd<N1>((T*) po.data());
		for (int k1 = 0; k1 < N1o2; k1++) {
			const int k = N2 * k1 + k2;
			X[k] = pe[k1].real() + po[k1].real();
			X[N - k] = pe[k1].imag() + po[k1].imag();
			assert(N - k > k);
		}
		for (int k1 = N1o2; k1 < N1; k1++) {
			const int k = N2 * k1 + k2;
			X[N - k] = pe[k1 - N1o2].real() - po[k1 - N1o2].real();
			X[k] = -pe[k1 - N1o2].imag() + po[k1 - N1o2].imag();
			assert(N - k < k);
		}
	}
	if (N2 % 2 == 0) {
		std::array<T, N1o2> qe;
		std::array<T, N1o2> qo;
		for (int n1 = 0; n1 < N1o2; n1++) {
			qe[n1] = X[N2 * n1 + N2o2];
			qo[n1] = X[No2 + N2 * n1 + N2o2];
		}
		sfft_skew_odd<N1>(qo.data());
		for (int k1 = 0; k1 < N1o4; k1++) {
			X[N2 * k1 + N2o2] = qe[k1] + qo[k1];
			X[N - (N2 * k1 + N2o2)] = qe[N1o2 - k1 - 1] + qo[N1o2 - k1 - 1];
		}
		for (int k1 = N1o4; k1 < N1o2; k1++) {
			X[N2 * k1 + N2o2] = qe[N1o2 - k1 - 1] - qo[N1o2 - k1 - 1];
			X[N - (N2 * k1 + N2o2)] = -qe[k1] + qo[k1];
		}
	}
}

void fft_radix2(double* X, int N) {
	int j = 0;
	for (int i = 0; i < N - 1; i++) {
		if (i > j) {
			std::swap(X[i], X[j]);
		}
		int k = N / 2;
		while (k <= j) {
			j -= k;
			k >>= 1;
		}
		j += k;
	}
	fft_radix2<double>(X, N);
}

