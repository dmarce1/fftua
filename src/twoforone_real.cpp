#include "fftu.hpp"
#include "util.hpp"

void fft_twoforone_real(fft_simd4* X, int N) {
	complex<fft_simd4>* Z = (complex<fft_simd4>*) X;
	static std::vector<complex<fft_simd4>> Xe, Xo;
	const int No2 = N / 2;
	const int No4 = N / 4;
	const auto& W = twiddles(N);
	Xe.resize(No4 + 1);
	Xo.resize(No4 + 1);

	fft_scramble(Z, No2);
	fft(Z, No2);

	Xe[0].real() = Z[0].real();
	Xo[0].real() = Z[0].imag();
	Xe[0].imag() = Xo[0].imag() = 0.0;
	for (int k = 1; k < (No2 + 1) / 2; k++) {
		Xe[k] = 0.5 * (Z[k] + Z[No2 - k].conj());
		Xo[k] = 0.5 * (Z[No2 - k] - Z[k].conj());
		std::swap(Xo[k].real(), Xo[k].imag());
	}
	if (No2 % 2 == 0) {
		Xe[No2 / 2] = 0.5 * (Z[No4] + Z[No4].conj());
		Xo[No4] = 0.5 * (Z[No4] - Z[No4].conj());
		std::swap(Xo[No4].real(), Xo[No4].imag());
	}
	for (int k = 1; k <= No4; k++) {
		Xo[k] *= W[k];
	}
	X[0] = Xe[0].real() + Xo[0].real();
	for (int k = 1; k < (No2 + 1) / 2; k++) {
		complex<fft_simd4> xk;
		xk = Xe[k] + Xo[k];
		X[k] = xk.real();
		X[N - k] = xk.imag();
	}
	if (No2 % 2 == 0) {
		complex<fft_simd4> xk;
		xk = Xe[No4] + Xo[No4];
		X[No4] = xk.real();
		X[N - No4] = xk.imag();
	}
	for (int k = No4 + 1; k < No2; k++) {
		complex<fft_simd4> xk;
		xk = (Xe[No2 - k] - Xo[No2 - k]).conj();
		X[k] = xk.real();
		X[N - k] = xk.imag();
	}
	X[No2] = Xe[0].real() - Xo[0].real();
}

void fft_twoforone_real(double* X, int N) {
	complex<double>* Z = (complex<double>*) X;
	std::vector<complex<double>> Xe, Xo;
	const int No2 = N / 2;
	const auto& W = twiddles(N);
	Xe.resize(No2 / 2 + SIMD_SIZE);
	Xo.resize(No2 / 2 + SIMD_SIZE);

	fft(Z, No2);

	Xe[0] = 0.5 * (Z[0] + Z[0].conj());
	Xo[0] = 0.5 * (Z[0] - Z[0].conj());
	std::swap(Xo[0].real(), Xo[0].imag());
	Xe[1] = 0.5 * (Z[1] + Z[No2 - 1].conj());
	Xo[1] = 0.5 * (Z[No2 - 1] - Z[1].conj());
	std::swap(Xo[1].real(), Xo[1].imag());
	int end = ((No2 / 2 + 1) % 2 == 0) ? (No2 / 2 + 1) : (No2 / 2);
	__m256d half = { 0.5, 0.5, 0.5, 0.5 };
	__m256d conj = { 1.0, -1.0, 1.0, -1.0 };
	for (int k = 2; k < end; k += 2) {
		__m256d zlo = *((__m256d *) (&Z[k]));
		__m256d zhi;
		zhi[0] = Z[No2 - k].real();
		zhi[1] = Z[No2 - k].imag();
		zhi[2] = Z[No2 - k - 1].real();
		zhi[3] = Z[No2 - k - 1].imag();
		__m256d& xe = *((__m256d*)(&Xe[k]));
		__m256d& xo = *((__m256d*)(&Xo[k]));
		xe = _mm256_mul_pd(half, _mm256_add_pd(zlo, _mm256_mul_pd(conj, zhi)));
		xo = _mm256_mul_pd(half, _mm256_sub_pd(zhi, _mm256_mul_pd(conj, zlo)));
		xo = _mm256_permute4x64_pd(xo, (1 << 0) | (0 << 2) | (3 << 4) | (2 << 6));
	}
	if (end == No2 / 2) {
		Xe[No2 / 2] = 0.5 * (Z[No2 / 2] + Z[No2 - No2 / 2].conj());
		Xo[No2 / 2] = 0.5 * (Z[No2 - No2 / 2] - Z[No2 / 2].conj());
		std::swap(Xo[No2 / 2].real(), Xo[No2 / 2].imag());
	}
	end = ((No2 / 2) % 2 == 1) ? No2 / 2 + 1 : No2 / 2;
	for (int k = 0; k < end; k += 2) {
		*((__m256d *) (&Xo[k])) = mul(*((__m256d *) (&Xo[k])), *((__m256d *) (&W[k])));
	}
	if (end == No2 / 2) {
		Xo[No2 / 2] *= W[No2 / 2];
	}
	complex<double> xk;
	X[0] = Xe[0].real() + Xo[0].real();
	for (int k = 1; k < (No2 + 1) / 2; k++) {
		xk = Xe[k] + Xo[k];
		X[k] = xk.real();
		X[N - k] = xk.imag();
	}
	if (No2 % 2 == 0) {
		xk = Xe[No2 / 2] + Xo[No2 / 2];
		X[No2 / 2] = xk.real();
		X[N - No2 / 2] = xk.imag();
	}
	for (int k = No2 / 2 + 1; k < No2; k++) {
		xk = (Xe[No2 - k] - Xo[No2 - k]).conj();
		X[k] = xk.real();
		X[N - k] = xk.imag();
	}
	X[N / 2] = Xe[0].real() - Xo[0].real();
}
