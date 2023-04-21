#include "fftu.hpp"
#include "util.hpp"

template<class T>
void fft_twoforone_real1(T* X, int N) {
	complex<T>* Z = (complex<T>*) X;
	static std::vector<complex<T>> Xe, Xo;
	const int No2 = N / 2;
	const auto& W = twiddles(N);
	Xe.resize(No2 / 2 + 1);
	Xo.resize(No2 / 2 + 1);

	fft(Z, No2);

	Xe[0] = 0.5 * (Z[0] + Z[0].conj());
	Xo[0] = 0.5 * (Z[0] - Z[0].conj());
	std::swap(Xo[0].real(), Xo[0].imag());
	for (int k = 1; k < (No2 + 1) / 2; k++) {
		Xe[k] = 0.5 * (Z[k] + Z[No2 - k].conj());
		Xo[k] = 0.5 * (Z[No2 - k] - Z[k].conj());
		std::swap(Xo[k].real(), Xo[k].imag());
	}
	if (No2 % 2 == 0) {
		Xe[No2 / 2] = 0.5 * (Z[No2 / 2] + Z[No2 / 2].conj());
		Xo[No2 / 2] = 0.5 * (Z[No2 / 2] - Z[No2 / 2].conj());
		std::swap(Xo[No2 / 2].real(), Xo[No2 / 2].imag());
	}
	X[0] = Xe[0].real() + Xo[0].real();
	for (int k = 1; k < (No2 + 1) / 2; k++) {
		complex<T> xk;
		xk = Xe[k] + Xo[k] * W[k];
		X[k] = xk.real();
		X[N - k] = xk.imag();
	}
	if (No2 % 2 == 0) {
		complex<T> xk;
		xk = Xe[No2 / 2] + Xo[No2 / 2] * W[No2 / 2];
		X[No2 / 2] = xk.real();
		X[N - No2 / 2] = xk.imag();
	}
	for (int k = No2 / 2 + 1; k < No2; k++) {
		complex<T> xk;
		xk = Xe[No2 - k].conj() + Xo[No2 - k].conj() * W[k];
		X[k] = xk.real();
		X[N - k] = xk.imag();
	}
	X[N / 2] = Xe[0].real() - Xo[0].real();
}

void fft_twoforone_real(fft_simd4* X, int N) {
	fft_twoforone_real1(X, N);
}

void fft_twoforone_real(double* X, int N) {
	fft_twoforone_real1(X, N);
}
