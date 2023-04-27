#include "sfft.hpp"
#include "util.hpp"

#include <array>
#include <cmath>
#include <cstring>
#include <numeric>
#include "fftu.hpp"
#include <stack>

void fft_split_real_indices(int R, int* I, int N) {
	if (N <= SFFT_NMAX) {
		return;
	}
	if (!(1 << ilogb(N) == N && R == 4)) {
		const int N1 = R;
		const int N1o2 = N1 / 2;
		const int N2 = N / N1;
		const int No2 = N / 2;
		std::vector<int> J;
		J.resize(N);
		for (int n = 0; n < No2; n++) {
			J[n] = I[2 * n];
		}
		for (int n1 = 0; n1 < N1o2; n1++) {
			for (int n2 = 0; n2 < N2; n2++) {
				J[No2 + N2 * n1 + n2] = I[N1 * n2 + 2 * n1 + 1];
			}
		}
		fft_indices_real(J.data(), No2);
		for (int n1 = 0; n1 < N1o2; n1++) {
			const int o = No2 + N2 * n1;
			fft_indices_real(J.data() + o, N2);
		}
		std::memcpy(I, J.data(), sizeof(int) * N);
	}
}

template<class T, int N1>
void fft_split_real(T* X, int N) {
	if (N <= SFFT_NMAX) {
		sfft_real(X, N);
		return;
	}
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

template<class T>
void fft_split_real_2pow(T* X, int N) {

	if (N <= 256) {
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

void fft_2pow(double* X, int N) {
	if (N <= 64) {
		sfft_real(X, N);
		return;
	}
	/*	const int No8 = N / 8;
	 const int No4 = N / 4;
	 const int No2 = N / 2;
	 int j = 0;
	 const int Nm1 = N - 1;
	 for (int i = 0; i < Nm1; i++) {
	 if (i < j) {
	 std::swap(X[i], X[j]);
	 }
	 int k = No2;
	 while (k <= j) {
	 j -= k;
	 k >>= 1;
	 }
	 j += k;
	 }
	 fft_split_real_2pow(X, N);*/
	constexpr int R = 4;
	constexpr int Ro2 = 2;
	fft_simd4* Z = (fft_simd4*) X;
	static std::vector<double> Y;
	Y.resize(N);
	const int No8 = N / 8;
	const int No4 = N / 4;
	const int No2 = N / 2;
	int j = 0;
	const int Nm1 = No4 - 1;
	for (int i = 0; i < Nm1; i++) {
		if (i < j) {
			std::swap(Z[i], Z[j]);
		}
		int k = No8;
		while (k <= j) {
			j -= k;
			k >>= 1;
		}
		j += k;
	}
	fft_split_real_2pow(Z, No4);
	for (int k2 = 0; k2 < No4; k2++) {
		for (int n1 = 0; n1 < R; n1++) {
			Y[No4 * n1 + k2] = X[R * k2 + n1];
		}
	}
	const auto& W = twiddles(N);
	const auto& Wv = vector_twiddles(No4, R);
	std::array<double, R> q;
	q[0] = Y[0];
	q[1] = Y[No4];
	q[2] = Y[No2];
	q[3] = Y[No2 + No4];
	sfft_real<R>(q.data());
	X[0] = q[0];
	X[No4] = q[1];
	X[No2] = q[2];
	X[No2 + No4] = q[3];
	for (int n1 = 0; n1 < R; n1++) {
		q[n1] = Y[No4 * n1 + No8];
	}
	sfft_skew<R>(q.data());
	for (int k1 = 0; k1 < Ro2; k1++) {
		X[0 + (No4 * k1 + No8)] = q[k1];
		X[N - (No4 * k1 + No8)] = q[R - k1 - 1];
	}
	for (int k2 = 1; k2 < std::min(No8, SIMD_SIZE); k2++) {
		std::array<complex<double>, R> p;
		for (int n1 = 0; n1 < R; n1++) {
			p[n1].real() = Y[No4 * n1 + k2];
			p[n1].imag() = Y[No4 * n1 - k2 + No4];
			p[n1] *= W[n1 * k2];
		}
		sfft_complex<R>((double*) p.data());
		for (int k1 = 0; k1 < Ro2; k1++) {
			const int k = No4 * k1 + k2;
			X[k] = p[k1].real();
			X[N - k] = p[k1].imag();
		}
		for (int k1 = Ro2; k1 < R; k1++) {
			const int k = No4 * k1 + k2;
			X[N - k] = p[k1].real();
			X[k] = -p[k1].imag();
		}
	}
	alignas(32)
	std::array<complex<fft_simd4>, R> p;
	for (int k2 = SIMD_SIZE; k2 < No8; k2 += SIMD_SIZE) {
		for (int n1 = 0; n1 < R; n1++) {
			p[n1].real() = *((fft_simd4*) (&Y[No4 * n1 + k2]));
			p[n1].imag().load(&Y[No4 * n1 + No4 - k2 - SIMD_SIZE + 1]);
			p[n1].imag() = p[n1].imag().reverse();
		}
		for (int n1 = 1; n1 < R; n1++) {
			p[n1] *= Wv[k2 / SIMD_SIZE][n1];
		}
		sfft_complex<R>((fft_simd4*) p.data());
		for (int k1 = 0; k1 < Ro2; k1++) {
			const int k = No4 * k1 + k2;
			*((fft_simd4*) &X[k]) = p[k1].real();
			for (int i = 0; i < SIMD_SIZE; i++) {
				X[N - k - i] = p[k1].imag()[i];
			}
		}
		for (int k1 = Ro2; k1 < R; k1++) {
			const int k = No4 * k1 + k2;
			for (int i = 0; i < SIMD_SIZE; i++) {
				X[N - k - i] = p[k1].real()[i];
			}
			*((fft_simd4*) &X[k]) = -p[k1].imag();
		}
	}
}

template<class T>
void fft_split1_real(int N1, T* X, int N) {
	switch (N1) {
	case 4:
		if (1 << ilogb(N) == N) {
			return fft_split_real_2pow<T>(X, N);
		} else {
			return fft_split_real<T, 4>(X, N);
		}
	case 8:
		return fft_split_real<T, 8>(X, N);
	case 12:
		return fft_split_real<T, 12>(X, N);
	case 16:
		return fft_split_real<T, 16>(X, N);
	case 20:
		return fft_split_real<T, 20>(X, N);
	case 24:
		return fft_split_real<T, 24>(X, N);
	case 28:
		return fft_split_real<T, 28>(X, N);
	case 32:
		return fft_split_real<T, 32>(X, N);
	case 36:
		return fft_split_real<T, 36>(X, N);
	case 40:
		return fft_split_real<T, 40>(X, N);
	case 44:
		return fft_split_real<T, 44>(X, N);
	case 48:
		return fft_split_real<T, 48>(X, N);
	case 52:
		return fft_split_real<T, 52>(X, N);
	case 56:
		return fft_split_real<T, 56>(X, N);
	case 60:
		return fft_split_real<T, 60>(X, N);
	case 64:
		return fft_split_real<T, 64>(X, N);
	}
}

void fft_split_real(int N1, double* X, int N) {
	fft_split1_real(N1, X, N);
}

void fft_split_real(int N1, fft_simd4* X, int N) {
	fft_split1_real(N1, X, N);
}
