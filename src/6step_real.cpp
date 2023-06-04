#include "fftu.hpp"
#include "util.hpp"
#include <memory>

int bit_reverse(int i, int N) {
	const auto w = ilogb(N);
	int j = 0;
	for (int k = 0; k < w; k++) {
		j <<= 1;
		j |= i & 1;
		i >>= 1;
	}
	return j;
}

const std::vector<int>& bit_reversals(int N) {
	using entry_type = std::shared_ptr<std::vector<int>>;
	static std::unordered_map<int, entry_type> cache;
	auto iter = cache.find(N);
	if (iter != cache.end()) {
		return *(iter->second);
	} else {
		std::vector<int> b(N);
		for (int n = 0; n < N; n++) {
			b[n] = bit_reverse(n, N);
		}
		cache[N] = std::make_shared<std::vector<int>>(std::move(b));
		return *(cache[N]);
	}
}

#include <cstring>

template<class T>
void scramble(T* x, int N) {
	static std::vector<T> buffer;
	const int M = lround(sqrt(N));
	buffer.resize(M);
	int b = 4;
	const auto& BR = bit_reversals(M);
	for (int n0 = 0; n0 < M; n0 += b) {
		int m0 = n0;
		const int n1 = std::min(n0 + b, M);
		for (int n = n0; n < n1; n++) {
			const int m1 = std::min(m0 + b, M);
			for (int m = n + 1; m < m1; m++) {
				std::swap(x[M * m + n], x[M * n + m]);
			}
		}
		const auto mask = _mm_set_epi32(12, 8, 4, 0);
		for (int m0 = n0 + b; m0 < M; m0 += b) {
			__m256d A[4];
			__m256d B[4];
			A[0] = _mm256_load_pd(x + M * (n0 + 0) + m0);
			A[1] = _mm256_load_pd(x + M * (n0 + 1) + m0);
			A[2] = _mm256_load_pd(x + M * (n0 + 2) + m0);
			A[3] = _mm256_load_pd(x + M * (n0 + 3) + m0);
			B[0] = _mm256_i32gather_pd((double* ) A + 0, mask, 8);
			B[1] = _mm256_i32gather_pd((double* ) A + 1, mask, 8);
			B[2] = _mm256_i32gather_pd((double* ) A + 2, mask, 8);
			B[3] = _mm256_i32gather_pd((double* ) A + 3, mask, 8);
			A[0] = _mm256_load_pd(x + M * (m0 + 0) + n0);
			A[1] = _mm256_load_pd(x + M * (m0 + 1) + n0);
			A[2] = _mm256_load_pd(x + M * (m0 + 2) + n0);
			A[3] = _mm256_load_pd(x + M * (m0 + 3) + n0);
			_mm256_store_pd(x + M * (m0 + 0) + n0, B[0]);
			_mm256_store_pd(x + M * (m0 + 1) + n0, B[1]);
			_mm256_store_pd(x + M * (m0 + 2) + n0, B[2]);
			_mm256_store_pd(x + M * (m0 + 3) + n0, B[3]);
			B[0] = _mm256_i32gather_pd((double* ) A + 0, mask, 8);
			B[1] = _mm256_i32gather_pd((double* ) A + 1, mask, 8);
			B[2] = _mm256_i32gather_pd((double* ) A + 2, mask, 8);
			B[3] = _mm256_i32gather_pd((double* ) A + 3, mask, 8);
			_mm256_store_pd(x + M * (n0 + 0) + m0, B[0]);
			_mm256_store_pd(x + M * (n0 + 1) + m0, B[1]);
			_mm256_store_pd(x + M * (n0 + 2) + m0, B[2]);
			_mm256_store_pd(x + M * (n0 + 3) + m0, B[3]);
		}
	}
	for (int i = 0; i < M; i++) {
		int j = BR[i];
		if (i < j) {
			std::memcpy(buffer.data(), x + M * i, sizeof(T) * M);
			std::memcpy(x + M * i, x + M * j, sizeof(T) * M);
			std::memcpy(x + M * j, buffer.data(), sizeof(T) * M);
		}
	}
}

void fft_batch_real1(double* x, int M, int N) {
	constexpr int N1 = 4;
	const int N2 = N / N1;
	std::array<complex<fft_simd4>, SFFT_NMAX> z0;
	const auto& W = twiddles(N);
	__m256d R0, R1, R2, R3;
	__m256d S0, S1, S2, S3;
	__m256d T0, T1, T2, T3;
	__m256d U0, U1, U2, U3;
	const __m256d C0 = _mm256_set1_pd(M_SQRT1_2);
	const __m256d C1 = _mm256_set1_pd(-M_SQRT1_2);
	double *xi0, *xi1, *xi2, *xi3, *xi4, *xi5, *xi6, *xi7;
	const int D = M * N2;
	if (N == 2) {
		for (int m = 0; m < M; m += SIMD_SIZE) {
			int i0 = 0;
			i0 = m + M * i0;
			xi0 = x + i0;
			xi1 = xi0 + D;
			R0 = _mm256_load_pd(xi0);
			R1 = _mm256_load_pd(xi1);
			R0 = _mm256_add_pd(R0, R1);
			R1 = _mm256_sub_pd(R0, R1);
			_mm256_store_pd(xi0, R0);
			_mm256_store_pd(xi1, R1);
		}
		return;
	}
	if (N <= 1) {
		return;
	}
	for (int n1 = 0; n1 < N1; n1++) {
		fft_batch_real1(x + n1 * N2 * M, M, N2);
	}
	for (int m = 0; m < M; m += SIMD_SIZE) {
		int i0 = 0;
		i0 = m + M * i0;
		xi0 = x + i0;
		xi1 = xi0 + D;
		xi2 = xi1 + D;
		xi3 = xi2 + D;
		R0 = _mm256_load_pd(xi0);
		R1 = _mm256_load_pd(xi2);
		R2 = _mm256_load_pd(xi1);
		R3 = _mm256_load_pd(xi3);
		S0 = _mm256_add_pd(R0, R2);
		S1 = _mm256_sub_pd(R0, R2);
		S2 = _mm256_add_pd(R1, R3);
		S3 = _mm256_sub_pd(R3, R1);
		R0 = _mm256_add_pd(S0, S2);
		R1 = _mm256_sub_pd(S0, S2);
		_mm256_store_pd(xi0, R0);
		_mm256_store_pd(xi1, S1);
		_mm256_store_pd(xi2, R1);
		_mm256_store_pd(xi3, S3);
	}
	if (N2 > 1) {
		for (int m = 0; m < M; m += SIMD_SIZE) {
			int i0 = N2 / 2;
			i0 = m + M * i0;
			xi0 = x + i0;
			xi1 = xi0 + D;
			xi2 = xi1 + D;
			xi3 = xi2 + D;
			S0 = _mm256_load_pd(xi0);
			S1 = _mm256_load_pd(xi2);
			S3 = _mm256_load_pd(xi1);
			S2 = _mm256_load_pd(xi3);
			R1 = _mm256_sub_pd(S1, S2);
			R2 = _mm256_add_pd(S1, S2);
			S1 = _mm256_mul_pd(R1, C0);
			S2 = _mm256_mul_pd(R2, C1);
			T0 = _mm256_add_pd(S0, S1);
			T1 = _mm256_sub_pd(S2, S3);
			T2 = _mm256_sub_pd(S0, S1);
			T3 = _mm256_add_pd(S3, S2);
			_mm256_store_pd(xi0, T0);
			_mm256_store_pd(xi1, T2);
			_mm256_store_pd(xi2, T3);
			_mm256_store_pd(xi3, T1);
		}
	}
	for (int k2 = 1; k2 < N2 / 2; k2++) {
		const int j = k2;
		const auto w1 = W[1 * j];
		const auto w2 = W[2 * j];
		const auto w3 = W[3 * j];
		const __m256d cos1 = _mm256_set1_pd(w1.real());
		const __m256d cos2 = _mm256_set1_pd(w2.real());
		const __m256d cos3 = _mm256_set1_pd(w3.real());
		const __m256d sin1 = _mm256_set1_pd(w1.imag());
		const __m256d sin2 = _mm256_set1_pd(w2.imag());
		const __m256d sin3 = _mm256_set1_pd(w3.imag());
		for (int m = 0; m < M; m += SIMD_SIZE) {
			int i0 = k2;
			int i1 = N2 - k2;
			i0 = m + M * i0;
			i1 = m + M * i1;
			xi0 = x + i0;
			xi1 = x + i1;
			xi2 = xi0 + D;
			xi3 = xi1 + D;
			xi4 = xi2 + D;
			xi5 = xi3 + D;
			xi6 = xi4 + D;
			xi7 = xi5 + D;
			R0 = _mm256_load_pd(xi0);
			S0 = _mm256_load_pd(xi1);
			R1 = _mm256_load_pd(xi4);
			S1 = _mm256_load_pd(xi5);
			R2 = _mm256_load_pd(xi2);
			S2 = _mm256_load_pd(xi3);
			R3 = _mm256_load_pd(xi6);
			S3 = _mm256_load_pd(xi7);
			T0 = _mm256_mul_pd(S1, sin1);
			T1 = _mm256_mul_pd(S1, cos1);
			S1 = _mm256_fmadd_pd(sin1, R1, T1);
			R1 = _mm256_fmsub_pd(cos1, R1, T0);
			T0 = _mm256_mul_pd(S2, sin2);
			T1 = _mm256_mul_pd(S2, cos2);
			S2 = _mm256_fmadd_pd(sin2, R2, T1);
			R2 = _mm256_fmsub_pd(cos2, R2, T0);
			T0 = _mm256_mul_pd(S3, sin3);
			T1 = _mm256_mul_pd(S3, cos3);
			S3 = _mm256_fmadd_pd(sin3, R3, T1);
			R3 = _mm256_fmsub_pd(cos3, R3, T0);
			U0 = _mm256_add_pd(R0, R2);
			T0 = _mm256_add_pd(S0, S2);
			U2 = _mm256_sub_pd(R0, R2);
			T2 = _mm256_sub_pd(S0, S2);
			U1 = _mm256_add_pd(R1, R3);
			T1 = _mm256_add_pd(S1, S3);
			U3 = _mm256_sub_pd(R3, R1);
			T3 = _mm256_sub_pd(S1, S3);
			R0 = _mm256_add_pd(U0, U1);
			S0 = _mm256_add_pd(T0, T1);
			R1 = _mm256_add_pd(U2, T3);
			S1 = _mm256_add_pd(T2, U3);
			R2 = _mm256_sub_pd(U0, U1);
			S2 = _mm256_sub_pd(T1, T0);
			R3 = _mm256_sub_pd(U2, T3);
			S3 = _mm256_sub_pd(U3, T2);
			_mm256_store_pd(xi0, R0);
			_mm256_store_pd(xi7, S0);
			_mm256_store_pd(xi2, R1);
			_mm256_store_pd(xi5, S1);
			_mm256_store_pd(xi3, R2);
			_mm256_store_pd(xi4, S2);
			_mm256_store_pd(xi1, R3);
			_mm256_store_pd(xi6, S3);
		}
	}
}

void fft_batch_real(double* x, int L, int N) {
	int N1 = (ilogb(N) % 2 == 0) ? 4 : 8;
	int N2 = 1;
	int NHI = N / (N1 * N2);
	const auto& W = twiddles(N);
	__m256d R0, R1, R2, R3, R4, R5, R6, R7;
	__m256d S0, S1, S2, S3, S6;
	__m256d T0, T1, T2, T3;
	__m256d U0, U1, U2, U3, U4;
	double *xi0, *xi1, *xi2, *xi3, *xi4, *xi5, *xi6, *xi7;
	const __m256d s1 = _mm256_set1_pd(-M_SQRT1_2);
	const __m256d c1 = _mm256_set1_pd(+M_SQRT1_2);
	const __m256d C0 = _mm256_set1_pd(M_SQRT1_2);
	const __m256d C1 = _mm256_set1_pd(-M_SQRT1_2);
	if (N1 == 8) {
		int D = L * N2;
		for (int ihi = 0; ihi < NHI; ihi++) {
			for (int l = 0; l < L; l += SIMD_SIZE) {
				int i0 = N2 * N1 * ihi;
				i0 = l + L * i0;
				xi0 = x + i0;
				xi1 = xi0 + D;
				xi2 = xi1 + D;
				xi3 = xi2 + D;
				xi4 = xi3 + D;
				xi5 = xi4 + D;
				xi6 = xi5 + D;
				xi7 = xi6 + D;
				R0 = _mm256_load_pd(xi0);
				R1 = _mm256_load_pd(xi4);
				R2 = _mm256_load_pd(xi2);
				R3 = _mm256_load_pd(xi6);
				R4 = _mm256_load_pd(xi1);
				R5 = _mm256_load_pd(xi5);
				R6 = _mm256_load_pd(xi3);
				R7 = _mm256_load_pd(xi7);
				U0 = _mm256_add_pd(R0, R4);
				U2 = _mm256_sub_pd(R0, R4);
				U1 = _mm256_add_pd(R2, R6);
				U3 = _mm256_sub_pd(R2, R6);
				R0 = _mm256_add_pd(U0, U1);
				R6 = R2 = U2;
				S6 = S2 = U3;
				R4 = _mm256_sub_pd(U0, U1);
				U0 = _mm256_add_pd(R1, R5);
				U2 = _mm256_sub_pd(R1, R5);
				U1 = _mm256_add_pd(R3, R7);
				U3 = _mm256_sub_pd(R3, R7);
				R1 = _mm256_add_pd(U0, U1);
				R5 = _mm256_sub_pd(U0, U1);
				T0 = _mm256_mul_pd(U3, c1);
				T1 = _mm256_mul_pd(U3, s1);
				S3 = _mm256_fmadd_pd(s1, U2, T1);
				R3 = _mm256_fmsub_pd(c1, U2, T0);
				U0 = _mm256_add_pd(R0, R1);
				U4 = _mm256_sub_pd(R0, R1);
				U1 = _mm256_add_pd(R2, R3);
				T1 = _mm256_sub_pd(S3, S2);
				U3 = _mm256_sub_pd(R6, R3);
				T3 = _mm256_add_pd(S6, S3);
				_mm256_store_pd(xi0, U0);
				_mm256_store_pd(xi1, U1);
				_mm256_store_pd(xi2, R4);
				_mm256_store_pd(xi3, U3);
				_mm256_store_pd(xi4, U4);
				_mm256_store_pd(xi5, T3);
				_mm256_store_pd(xi6, -R5);
				_mm256_store_pd(xi7, T1);

			}
		}
		N2 *= N1;
		N1 = 4;
		NHI = N / (N1 * N2);
	}
	while (N2 < N) {
		int D = L * N2;
		for (int ihi = 0; ihi < NHI; ihi++) {
			for (int l = 0; l < L; l += SIMD_SIZE) {
				int i0 = N2 * N1 * ihi;
				i0 = l + L * i0;
				xi0 = x + i0;
				xi1 = xi0 + D;
				xi2 = xi1 + D;
				xi3 = xi2 + D;
				R0 = _mm256_load_pd(xi0);
				R1 = _mm256_load_pd(xi2);
				R2 = _mm256_load_pd(xi1);
				R3 = _mm256_load_pd(xi3);
				S0 = _mm256_add_pd(R0, R2);
				S1 = _mm256_sub_pd(R0, R2);
				S2 = _mm256_add_pd(R1, R3);
				S3 = _mm256_sub_pd(R3, R1);
				R0 = _mm256_add_pd(S0, S2);
				R1 = _mm256_sub_pd(S0, S2);
				_mm256_store_pd(xi0, R0);
				_mm256_store_pd(xi1, S1);
				_mm256_store_pd(xi2, R1);
				_mm256_store_pd(xi3, S3);
			}
		}
		if (N2 > 1) {
			for (int ihi = 0; ihi < NHI; ihi++) {
				for (int l = 0; l < L; l += SIMD_SIZE) {
					int i0 = N2 / 2 + N2 * N1 * ihi;
					i0 = l + L * i0;
					xi0 = x + i0;
					xi1 = xi0 + D;
					xi2 = xi1 + D;
					xi3 = xi2 + D;
					S0 = _mm256_load_pd(xi0);
					S1 = _mm256_load_pd(xi2);
					S3 = _mm256_load_pd(xi1);
					S2 = _mm256_load_pd(xi3);
					R1 = _mm256_sub_pd(S1, S2);
					R2 = _mm256_add_pd(S1, S2);
					S1 = _mm256_mul_pd(R1, C0);
					S2 = _mm256_mul_pd(R2, C1);
					T0 = _mm256_add_pd(S0, S1);
					T1 = _mm256_sub_pd(S2, S3);
					T2 = _mm256_sub_pd(S0, S1);
					T3 = _mm256_add_pd(S3, S2);
					_mm256_store_pd(xi0, T0);
					_mm256_store_pd(xi1, T2);
					_mm256_store_pd(xi2, T3);
					_mm256_store_pd(xi3, T1);
				}
			}
		}
		for (int k2 = 1; k2 < N2 / 2; k2++) {
			const int j = NHI * k2;
			const auto w1 = W[1 * j];
			const auto w2 = W[2 * j];
			const auto w3 = W[3 * j];
			const __m256d cos1 = _mm256_set1_pd(w1.real());
			const __m256d cos2 = _mm256_set1_pd(w2.real());
			const __m256d cos3 = _mm256_set1_pd(w3.real());
			const __m256d sin1 = _mm256_set1_pd(w1.imag());
			const __m256d sin2 = _mm256_set1_pd(w2.imag());
			const __m256d sin3 = _mm256_set1_pd(w3.imag());
			for (int ihi = 0; ihi < NHI; ihi++) {
				for (int l = 0; l < L; l += SIMD_SIZE) {
					int i0 = k2 + N2 * N1 * ihi;
					int i1 = N2 - k2 + N2 * N1 * ihi;
					i0 = l + L * i0;
					i1 = l + L * i1;
					xi0 = x + i0;
					xi1 = x + i1;
					xi2 = xi0 + D;
					xi3 = xi1 + D;
					xi4 = xi2 + D;
					xi5 = xi3 + D;
					xi6 = xi4 + D;
					xi7 = xi5 + D;
					R0 = _mm256_load_pd(xi0);
					S0 = _mm256_load_pd(xi1);
					R1 = _mm256_load_pd(xi4);
					S1 = _mm256_load_pd(xi5);
					R2 = _mm256_load_pd(xi2);
					S2 = _mm256_load_pd(xi3);
					R3 = _mm256_load_pd(xi6);
					S3 = _mm256_load_pd(xi7);
					T0 = _mm256_mul_pd(S1, sin1);
					T1 = _mm256_mul_pd(S1, cos1);
					S1 = _mm256_fmadd_pd(sin1, R1, T1);
					R1 = _mm256_fmsub_pd(cos1, R1, T0);
					T0 = _mm256_mul_pd(S2, sin2);
					T1 = _mm256_mul_pd(S2, cos2);
					S2 = _mm256_fmadd_pd(sin2, R2, T1);
					R2 = _mm256_fmsub_pd(cos2, R2, T0);
					T0 = _mm256_mul_pd(S3, sin3);
					T1 = _mm256_mul_pd(S3, cos3);
					S3 = _mm256_fmadd_pd(sin3, R3, T1);
					R3 = _mm256_fmsub_pd(cos3, R3, T0);
					U0 = _mm256_add_pd(R0, R2);
					T0 = _mm256_add_pd(S0, S2);
					U2 = _mm256_sub_pd(R0, R2);
					T2 = _mm256_sub_pd(S0, S2);
					U1 = _mm256_add_pd(R1, R3);
					T1 = _mm256_add_pd(S1, S3);
					U3 = _mm256_sub_pd(R3, R1);
					T3 = _mm256_sub_pd(S1, S3);
					R0 = _mm256_add_pd(U0, U1);
					S0 = _mm256_add_pd(T0, T1);
					R1 = _mm256_add_pd(U2, T3);
					S1 = _mm256_add_pd(T2, U3);
					R2 = _mm256_sub_pd(U0, U1);
					S2 = _mm256_sub_pd(T1, T0);
					R3 = _mm256_sub_pd(U2, T3);
					S3 = _mm256_sub_pd(U3, T2);
					_mm256_store_pd(xi0, R0);
					_mm256_store_pd(xi7, S0);
					_mm256_store_pd(xi2, R1);
					_mm256_store_pd(xi5, S1);
					_mm256_store_pd(xi3, R2);
					_mm256_store_pd(xi4, S2);
					_mm256_store_pd(xi1, R3);
					_mm256_store_pd(xi6, S3);
				}
			}
		}
		N2 *= N1;
		N1 = 4;
		NHI = N / (N1 * N2);
	}
}

#include <cstring>

extern "C" {
void fft_iter_real(double*, double*, int N, int M, __m256d dummy = __m256d());
}

template<class T>
void fft_6_real(T* x, int N) {
	/*{
		const auto& BR = bit_reversals(N);
		int M = 4;
		std::vector<double> xx(N * M);
		for (int n = 0; n < M; n++) {
			for (int m = 0; m < N; m++) {
				xx[n + M * m] = x[BR[m]];
			}
		}
		const auto& W3 = twiddles(N);
		fft_iter_real(xx.data(), (double*) W3.data(), N, M);
		for (int n = 0; n < M; n++) {
			for (int m = 0; m < N; m++) {
				x[m] = xx[n + M * m];
			}
		}
		return;
	}*/
	static std::vector<T> buffer;
	const int M = lround(sqrt(N));
	assert(M * M == N);
	const auto& W0 = twiddles(N);
	const auto& W1 = twiddles(2 * M);
	const auto& W2 = twiddles(M);
	const auto& BR = bit_reversals(M);
	const int No2 = N / 2;
	const int Mo2 = M / 2;
	T T1, T2, T3, T4;
	buffer.resize(M);

	for (int i = 0; i < M; i++) {
		int j = BR[i];
		if (i < j) {
			std::memcpy(buffer.data(), x + M * i, sizeof(T) * M);
			std::memcpy(x + M * i, x + M * j, sizeof(T) * M);
			std::memcpy(x + M * j, buffer.data(), sizeof(T) * M);
		}
	}
	//fft_batch_real(x, M, M);
	fft_iter_real(x, (double*) W2.data(), M, M);
	for (int k2 = 1; k2 < Mo2; k2++) {
		complex<double> wk2 = W0[k2];
		complex<__m256d> W, w0;
		W.real()[0] = 1.0;
		W.imag()[0] = 0.0;
		W.real()[1] = wk2.real();
		W.imag()[1] = wk2.imag();
		W.real()[2] = (wk2 * wk2).real();
		W.imag()[2] = (wk2 * wk2).imag();
		W.real()[3] = (wk2 * wk2 * wk2).real();
		W.imag()[3] = (wk2 * wk2 * wk2).imag();
		w0.real() = _mm256_set1_pd((wk2 * wk2 * wk2 * wk2).real());
		w0.imag() = _mm256_set1_pd((wk2 * wk2 * wk2 * wk2).imag());
		for (int n1 = 0; n1 < M; n1 += SIMD_SIZE) {
			const double* wptr = (double*) (&W);
			double* rptr = &x[M * k2 + n1];
			double* iptr = &x[M * (M - k2) + n1];
			const __m256d C = _mm256_load_pd(wptr + 0);
			const __m256d S = _mm256_load_pd(wptr + SIMD_SIZE);
			__m256d Re = _mm256_load_pd(rptr);
			__m256d Im = _mm256_load_pd(iptr);
			const __m256d T0 = _mm256_mul_pd(Re, S);
			const __m256d T1 = _mm256_mul_pd(Im, S);
			Im = _mm256_fmadd_pd(C, Im, T0);
			Re = _mm256_fmsub_pd(C, Re, T1);
			_mm256_store_pd(rptr, Re);
			_mm256_store_pd(iptr, Im);
			W *= w0;
		}
	}
	T1 = 0.0;
	T2 = 0.0;
	for (int n1 = 0; n1 < M; n1 += 2) {
		T1 += x[M * Mo2 + n1] * W1[n1].real();
		T2 += x[M * Mo2 + n1 + 1] * W1[n1 + 1].real();
	}
	for (int n1 = 0; n1 < M; n1++) {
		x[M * Mo2 + n1] *= -2.0 * W1[n1].imag();
	}

	scramble(x, N);

//	fft_batch_real(x, M, M);
	fft_iter_real(x, (double*) W2.data(), M, M);
	T3 = x[Mo2 + M * (M - 1)];
	x[Mo2 + M * (M - 1)] = -0.5 * x[Mo2];
	x[Mo2] = (T1 + T2);
	for (int k1 = 1; k1 < Mo2 - 1; k1++) {
		T4 = x[Mo2 + M * (M - k1 - 1)];
		x[Mo2 + M * (M - k1 - 1)] = x[Mo2 + M * (M - k1)] - x[Mo2 + M * k1];
		x[Mo2 + M * k1] = T3 + x[Mo2 + M * (k1 - 1)];
		T3 = T4;
	}
	x[Mo2 + No2] = 0.5 * x[No2 + Mo2];
	x[No2 - Mo2] = (T1 - T2);
	for (int k1 = 1; k1 < Mo2; k1++) {
		for (int k2 = 1; k2 < Mo2; k2++) {
			const int rx = M * k1 + k2;
			const int ry = N - M * k1 + k2;
			const int ix = M * k1 + M - k2;
			const int iy = N - M * k1 + M - k2;
			const auto rex = x[rx];
			const auto rey = x[ry];
			const auto imx = x[ix];
			const auto imy = x[iy];
			x[rx] = rex - imy;
			x[ix] = rey + imx;
			x[ry] = rex + imy;
			x[iy] = -rey + imx;
		}
	}

	for (int k1 = 0; k1 < Mo2; k1++) {
		for (int k2 = Mo2 + 1; k2 < M; k2++) {
			const int i = M * k1 + k2;
			const int j = N - M * k1 - k2;
			const int k = N - M * k1 - (M - k2);
			std::swap(x[i], x[j]);
			std::swap(x[j], x[k]);
			x[j] = -x[j];
		}
	}

}

void fft_6_real(double* x, int N) {
	fft_6_real<double>(x, N);
}
