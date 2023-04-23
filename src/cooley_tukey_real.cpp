#include "fftu.hpp"
#include "util.hpp"
#include <cstring>

void fft_cooley_tukey_indices_real(int N1, int* I, int N) {
	std::vector<int> J(N);
	const int N2 = N / N1;
	for (int n1 = 0; n1 < N1; n1++) {
		for (int n2 = 0; n2 < N2; n2++) {
			J[N2 * n1 + n2] = I[N1 * n2 + n1];
		}
	}
	std::memcpy(I, J.data(), sizeof(int) * N);
	for (int n1 = 0; n1 < N1; n1++) {
		fft_indices_real(&I[n1 * N2], N2);
	}
}

template<class T, int N1>
void fft_cooley_tukey_real(T* X, int N) {
	if (N <= SFFT_NMAX) {
		sfft_real(X, N);
		return;
	}

	constexpr int N1p1o2 = (N1 + 1) / 2;
	constexpr int N1o2 = N1 / 2;
	constexpr bool N1e = (N1 % 2 == 0);
	const int N2 = N / N1;
	const int No2 = N / 2;
	const int N2p1o2 = (N2 + 1) / 2;

	const auto& W = twiddles(N);

	for (int n1 = 0; n1 < N1; n1++) {
		fft_real(&X[n1 * N2], N2);
	}

	for (int k2 = 1; k2 < N2p1o2; k2++) {
		for (int n1 = 1; n1 < N1; n1++) {
			complex<T> z;
			const int ir = N2 * n1 + k2;
			const int ii = N2 * n1 - k2 + N2;
			z.real() = X[ir];
			z.imag() = X[ii];
			z *= W[n1 * k2];
			X[ir] = z.real();
			X[ii] = z.imag();
		}
	}

	std::array<T, N1> q;
	for (int n1 = 0; n1 < N1; n1++) {
		q[n1] = X[N2 * n1];
	}
	sfft_real<N1>(q.data());
	X[0] = q[0];
	for (int k1 = 1; k1 < N1p1o2; k1++) {
		X[0 + N2 * k1] = q[k1];
		X[N - N2 * k1] = q[N1 - k1];
	}
	if (N1e) {
		X[No2] = q[N1o2];
	}
	std::array<complex<T>, N1> p;
	for (int k2 = 1; k2 < N2p1o2; k2++) {
		for (int n1 = 0; n1 < N1; n1++) {
			p[n1].real() = X[N2 * n1 + k2];
			p[n1].imag() = X[N2 * n1 - k2 + N2];
		}
		sfft_complex<N1>((T*) p.data());
		for (int k1 = 0; k1 < N1p1o2; k1++) {
			const int k = N2 * k1 + k2;
			X[k] = p[k1].real();
			X[N - k] = p[k1].imag();
		}
		for (int k1 = N1p1o2; k1 < N1; k1++) {
			const int k = N2 * k1 + k2;
			X[N - k] = p[k1].real();
			X[k] = -p[k1].imag();
		}
	}
	if (N2 % 2 == 0) {
		std::array<T, N1> q;
		for (int n1 = 0; n1 < N1; n1++) {
			q[n1] = X[N2 * n1 + N2 / 2];
		}
		sfft_skew<N1>(q.data());
		for (int k1 = 0; k1 < N1p1o2; k1++) {
			X[0 + (N2 * k1 + N2 / 2)] = q[k1];
			X[N - (N2 * k1 + N2 / 2)] = q[N1 - k1 - 1];
		}
		if (!N1e) {
			X[No2] = q[N1o2];
		}
	}
}

template<class T>
void fft_cooley_tukey_real_big(int N1, T* X, int N) {
	if (N <= SFFT_NMAX) {
		sfft_real(X, N);
		return;
	}
	int N2 = N / N1;
	if( N2 % 2 == 0 ) {
		std::swap(N1, N2);
	}
	const int N1p1o2 = (N1 + 1) / 2;
	const int N1o2 = N1 / 2;
	const bool N1e = (N1 % 2 == 0);
	const int No2 = N / 2;
	const int N2p1o2 = (N2 + 1) / 2;

	const auto& W = twiddles(N);

	for (int n1 = 0; n1 < N1; n1++) {
		fft_real(&X[n1 * N2], N2);
	}

	for (int k2 = 1; k2 < N2p1o2; k2++) {
		for (int n1 = 1; n1 < N1; n1++) {
			complex<T> z;
			const int ir = N2 * n1 + k2;
			const int ii = N2 * n1 - k2 + N2;
			z.real() = X[ir];
			z.imag() = X[ii];
			z *= W[n1 * k2];
			X[ir] = z.real();
			X[ii] = z.imag();
		}
	}

	std::vector<T> q(N1);
	for (int n1 = 0; n1 < N1; n1++) {
		q[n1] = X[N2 * n1];
	}
	fft_scramble_real(q.data(), N1);
	fft_real(q.data(), N1);
	X[0] = q[0];
	for (int k1 = 1; k1 < N1p1o2; k1++) {
		X[0 + N2 * k1] = q[k1];
		X[N - N2 * k1] = q[N1 - k1];
	}
	if (N1e) {
		X[No2] = q[N1o2];
	}
	std::vector<complex<T>> p(N1);
	for (int k2 = 1; k2 < N2p1o2; k2++) {
		for (int n1 = 0; n1 < N1; n1++) {
			p[n1].real() = X[N2 * n1 + k2];
			p[n1].imag() = X[N2 * n1 - k2 + N2];
		}
		fft_scramble(p.data(), N1);
		fft(p.data(), N1);
		for (int k1 = 0; k1 < N1p1o2; k1++) {
			const int k = N2 * k1 + k2;
			X[k] = p[k1].real();
			X[N - k] = p[k1].imag();
		}
		for (int k1 = N1p1o2; k1 < N1; k1++) {
			const int k = N2 * k1 + k2;
			X[N - k] = p[k1].real();
			X[k] = -p[k1].imag();
		}
	}
}

template<class T>
void fft_cooley_tukey1_real(int N1, T* X, int N) {
	switch (N1) {
	case 2:
		return fft_cooley_tukey_real<T, 2>(X, N);
	case 3:
		return fft_cooley_tukey_real<T, 3>(X, N);
	case 4:
		return fft_cooley_tukey_real<T, 4>(X, N);
	case 5:
		return fft_cooley_tukey_real<T, 5>(X, N);
	case 6:
		return fft_cooley_tukey_real<T, 6>(X, N);
	case 7:
		return fft_cooley_tukey_real<T, 7>(X, N);
	case 8:
		return fft_cooley_tukey_real<T, 8>(X, N);
	case 9:
		return fft_cooley_tukey_real<T, 9>(X, N);
	case 10:
		return fft_cooley_tukey_real<T, 10>(X, N);
	case 11:
		return fft_cooley_tukey_real<T, 11>(X, N);
	case 12:
		return fft_cooley_tukey_real<T, 12>(X, N);
	case 13:
		return fft_cooley_tukey_real<T, 13>(X, N);
	case 14:
		return fft_cooley_tukey_real<T, 14>(X, N);
	case 15:
		return fft_cooley_tukey_real<T, 15>(X, N);
	case 16:
		return fft_cooley_tukey_real<T, 16>(X, N);
	case 17:
		return fft_cooley_tukey_real<T, 17>(X, N);
	case 18:
		return fft_cooley_tukey_real<T, 18>(X, N);
	case 19:
		return fft_cooley_tukey_real<T, 19>(X, N);
	case 20:
		return fft_cooley_tukey_real<T, 20>(X, N);
	case 21:
		return fft_cooley_tukey_real<T, 21>(X, N);
	case 22:
		return fft_cooley_tukey_real<T, 22>(X, N);
	case 23:
		return fft_cooley_tukey_real<T, 23>(X, N);
	case 24:
		return fft_cooley_tukey_real<T, 24>(X, N);
	case 25:
		return fft_cooley_tukey_real<T, 25>(X, N);
	case 26:
		return fft_cooley_tukey_real<T, 26>(X, N);
	case 27:
		return fft_cooley_tukey_real<T, 27>(X, N);
	case 28:
		return fft_cooley_tukey_real<T, 28>(X, N);
	case 29:
		return fft_cooley_tukey_real<T, 29>(X, N);
	case 30:
		return fft_cooley_tukey_real<T, 30>(X, N);
	case 31:
		return fft_cooley_tukey_real<T, 31>(X, N);
	case 32:
		return fft_cooley_tukey_real<T, 32>(X, N);
	case 33:
		return fft_cooley_tukey_real<T, 33>(X, N);
	case 34:
		return fft_cooley_tukey_real<T, 34>(X, N);
	case 35:
		return fft_cooley_tukey_real<T, 35>(X, N);
	case 36:
		return fft_cooley_tukey_real<T, 36>(X, N);
	case 37:
		return fft_cooley_tukey_real<T, 37>(X, N);
	case 38:
		return fft_cooley_tukey_real<T, 38>(X, N);
	case 39:
		return fft_cooley_tukey_real<T, 39>(X, N);
	case 40:
		return fft_cooley_tukey_real<T, 40>(X, N);
	case 41:
		return fft_cooley_tukey_real<T, 41>(X, N);
	case 42:
		return fft_cooley_tukey_real<T, 42>(X, N);
	case 43:
		return fft_cooley_tukey_real<T, 43>(X, N);
	case 44:
		return fft_cooley_tukey_real<T, 44>(X, N);
	case 45:
		return fft_cooley_tukey_real<T, 45>(X, N);
	case 46:
		return fft_cooley_tukey_real<T, 46>(X, N);
	case 47:
		return fft_cooley_tukey_real<T, 47>(X, N);
	case 48:
		return fft_cooley_tukey_real<T, 48>(X, N);
	case 49:
		return fft_cooley_tukey_real<T, 49>(X, N);
	case 50:
		return fft_cooley_tukey_real<T, 50>(X, N);
	case 51:
		return fft_cooley_tukey_real<T, 51>(X, N);
	case 52:
		return fft_cooley_tukey_real<T, 52>(X, N);
	case 53:
		return fft_cooley_tukey_real<T, 53>(X, N);
	case 54:
		return fft_cooley_tukey_real<T, 54>(X, N);
	case 55:
		return fft_cooley_tukey_real<T, 55>(X, N);
	case 56:
		return fft_cooley_tukey_real<T, 56>(X, N);
	case 57:
		return fft_cooley_tukey_real<T, 57>(X, N);
	case 58:
		return fft_cooley_tukey_real<T, 58>(X, N);
	case 59:
		return fft_cooley_tukey_real<T, 59>(X, N);
	case 60:
		return fft_cooley_tukey_real<T, 60>(X, N);
	case 61:
		return fft_cooley_tukey_real<T, 61>(X, N);
	case 62:
		return fft_cooley_tukey_real<T, 62>(X, N);
	case 63:
		return fft_cooley_tukey_real<T, 63>(X, N);
	case 64:
		return fft_cooley_tukey_real<T, 64>(X, N);
	default:
		return fft_cooley_tukey_real_big<T>(N1, X, N);
	}
}

void fft_cooley_tukey_real(int N1, double* X, int N) {
	fft_cooley_tukey1_real(N1, X, N);
}

void fft_cooley_tukey_real(int N1, fft_simd4* X, int N) {
	fft_cooley_tukey1_real(N1, X, N);
}

