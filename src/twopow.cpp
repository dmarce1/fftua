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
	if (N <= SFFT_NMAX) {
		scramble(X, N);
		sfft_real(X, N);
		return;
	}
	constexpr int N1 = 4;
	const int N2 = N / N1;
	const int N22 = N2 / 2;
	static std::vector<int> D;
	if (!D.size()) {
		D.resize(N1);
		std::iota(D.begin(), D.end(), 0);
		scramble(D.data(), N1);
	}
	for (int n1 = 0; n1 < N1; n1++) {
		fft2(&X[n1 * N2], N2);
	}
	const auto& C = cos_twiddles(N1, N2);
	const auto& S = sin_twiddles(N1, N2);
	T z[2 * N1];
	T* q = (T*) z;
	complex<T>* p = (complex<T>*) z;
	for (int n1 = 0; n1 < N1; n1++) {
		q[n1] = X[N2 * D[n1]];
	}
	sfft_real<N1>(q);
	X[0] = q[0];
	for (int k1 = 1; k1 < N1 / 2; k1++) {
		X[0 + N2 * k1] = q[k1];
		X[N - N2 * k1] = q[N1 - k1];
	}
	X[N / 2] = q[N1 / 2];
	if (N2 / 2) {
		for (int n1 = 0; n1 < N1; n1++) {
			q[n1] = X[N2 * D[n1] + N2 / 2];
		}
		sfft_skew<N1>(q);
		for (int k1 = 0; k1 < N1 / 2; k1++) {
			X[0 + (N2 * k1 + N2 / 2)] = q[k1];
			X[N - (N2 * k1 + N2 / 2)] = q[N1 - k1 - 1];
		}
	}
	for (int n1 = 1; n1 < N1; n1++) {
		const auto m1 = D[n1];
		const auto& c0 = C[n1];
		const auto& s0 = S[n1];
		const int N2m1 = N2 * m1;
		const int N2m1pN2 = N2m1 + N2;
		int k2 = 1;
		for (; k2 < N22; k2++) {
			auto& re = X[N2m1 + k2];
			auto& im = X[N2m1pN2 - k2];
			const auto re0 = re;
			re = re * c0[k2] - im * s0[k2];
			im = re0 * s0[k2] + im * c0[k2];
		}
	}
	for (int k2 = 1; k2 < N2 / 2; k2++) {
		for (int n1 = 0; n1 < N1; n1++) {
			p[n1].real() = X[N2 * D[n1] + k2];
			p[n1].imag() = X[N2 * D[n1] - k2 + N2];
		}
		sfft_complex<N1>((T*) p);
		for (int k1 = 0; k1 < N1 / 2; k1++) {
			const int k = N2 * k1 + k2;
			X[k] = p[k1].real();
			X[N - k] = p[k1].imag();
		}
		for (int k1 = N1 / 2; k1 < N1; k1++) {
			const int k = N2 * k1 + k2;
			X[N - k] = p[k1].real();
			X[k] = -p[k1].imag();
		}
	}
}

void fft_2pow(double* x, int N) {
	/*std::vector<int> I(N);
	 std::vector<int> J(N);
	 std::vector<bool> visited(N, false);
	 std::iota(I.begin(), I.end(), 0);
	 indices(I.data(), N);
	 J = I;
	 for (int n = 0; n < N; n++) {
	 I[J[n]] = n;
	 }
	 for (int n = 0; n < N; n++) {
	 printf("%i %i %i\n", n, I[n], reverse_bits_conj(n, N));
	 }
	 for (int n = 0; n < N; n++) {
	 if (!visited[n]) {
	 visited[n] = true;
	 auto tmp = x[n];
	 int m = reverse_bits_conj(n, N);
	 while (m != n) {
	 std::swap(tmp, x[m]);
	 visited[m] = true;
	 m = reverse_bits_conj(m, N);
	 }
	 std::swap(tmp, x[m]);
	 }
	 }*/
//	conjugate_scramble(x, N, N);
	scramble((fft_simd4*)x, N/4);
	fft3((fft_simd4*) x, N/4);
}

