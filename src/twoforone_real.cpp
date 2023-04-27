#include "fftu.hpp"
#include "util.hpp"
#include "types.hpp"

void fft_twoforone_real(fft_simd4* X, int N) {
	complex<fft_simd4>* Z = (complex<fft_simd4>*) X;
	const int No2 = N / 2;
	const int No4 = N / 4;
	const auto& W = twiddles(N);
	workspace<complex<fft_simd4>> ws;
	auto Xe = ws.create(No4 + 1);
	auto Xo = ws.create(No4 + 1);

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
	ws.destroy(std::move(Xo));
	ws.destroy(std::move(Xe));
}

void fft_twoforone_real(double* X, int N) {
	complex<double>* Z = (complex<double>*) X;
	const int No2 = N / 2;
	const auto& W = twiddles(N);
	workspace<complex<double>> ws;
	auto Xe = ws.create(No2 / 2 + SIMD_SIZE);
	auto Xo = ws.create(No2 / 2 + SIMD_SIZE);

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
	ws.destroy(std::move(Xo));
	ws.destroy(std::move(Xe));
}

void fft_fourfortwo_real(double* X, int N) {
	const int No2 = N / 2;
	const int No4 = N / 4;
	const auto& W = twiddles(N);
	workspace<complex<double>> ws;
	auto Xo = ws.create(No4);
	auto Xe = ws.create(No4);
	auto X0 = ws.create(No4 + 1);
	auto X1 = ws.create(No4 + 1);
	auto X2 = ws.create(No4 + 1);
	auto X3 = ws.create(No4 + 1);
	for (int n = 0; n < N / 4; n++) {
		Xe[n].real() = X[4 * n + 0];
		Xo[n].real() = X[4 * n + 1];
		Xe[n].imag() = X[4 * n + 2];
		Xo[n].imag() = X[4 * n + 3];
	}
	fft(Xe.data(), No4);
	fft(Xo.data(), No4);
	__m256d conj = { 1.0, -1.0, 1.0, -1.0 };
	__m256d half = { 0.5, 0.5, 0.5, 0.5 };
	constexpr int flip = (1 << 0) | (0 << 2) | (3 << 4) | (2 << 6) ;
	for (int k = 0; k < No4; k += (SIMD_SIZE / 2)) {
		auto xe0 = *((__m256d*) &Xe[k % No4]);
		auto xe1 = *((__m256d*) &Xe[(No4 - k) % No4]);
		auto xo0 = *((__m256d*) &Xo[k % No4]);
		auto xo1 = *((__m256d*) &Xo[(No4 - k) % No4]);
		xe1 = _mm256_permute4x64_pd(xe1, flip);
		xo1 = _mm256_permute4x64_pd(xo1, flip);
		auto z0 = _mm256_mul_pd(_mm256_add_pd(xe0, _mm256_mul_pd(xe1, conj)), half);
		auto z1 = _mm256_mul_pd(_mm256_add_pd(xo0, _mm256_mul_pd(xo1, conj)), half);
		auto z2 = _mm256_mul_pd(_mm256_add_pd(xe0, _mm256_mul_pd(xe1, conj)), half);
		auto z3 = _mm256_mul_pd(_mm256_add_pd(xo0, _mm256_mul_pd(xo1, conj)), half);
		z1 = _mm256_permute4x64_pd(z1, flip);
		z3 = _mm256_permute4x64_pd(z3, flip);
		z1 = _mm256_mul_pd(z1, conj);
		z2 = _mm256_mul_pd(z2, conj);
		__m256d tw1, tw2, tw3, us, ud, zs, zd, x0, x1, x2, x3;
		tw1[0] = W[k].real();
		tw1[1] = W[k].imag();
		tw1[2] = W[k].real();
		tw1[3] = W[k].imag();
		tw2[0] = W[2*k].real();
		tw2[1] = W[2*k].imag();
		tw2[2] = W[2*k].real();
		tw2[3] = W[2*k].imag();
		tw3[0] = W[3*k].real();
		tw3[1] = W[3*k].imag();
		tw3[2] = W[3*k].real();
		tw3[3] = W[3*k].imag();
		z1 = mul(z1, tw1);
		z2 = mul(z2, tw2);
		z3 = mul(z3, tw3);
		us = _mm256_add_pd(z0, z2);
		ud = _mm256_sub_pd(z0, z2);
		zs = _mm256_add_pd(z1, z3);
		zd = _mm256_sub_pd(z1, z3);
		z3 = _mm256_permute4x64_pd(z1, flip);
		x0 = _mm256_add_pd(us, zs);
		x1 = _mm256_sub_pd(ud, zd);
		x2 = _mm256_sub_pd(us, zs);
		x3 = _mm256_add_pd(ud, zd);
		X[N - (0 * N / 4 + k)] = x0[0];
		X[N - (1 * N / 4 + k)] = x1[0];
		X[2 * N / 4 + k] = -x2[0];
		X[3 * N / 4 + k] = -x3[0];
		X[0 * N / 4 + k] = x0[1];
		X[1 * N / 4 + k] = x1[1];
		X[N - (2 * N / 4 + k)] = x2[1];
		X[N - (3 * N / 4 + k)] = x3[1];
		X[N - (0 * N / 4 + k + 1)] = x0[2];
		X[N - (1 * N / 4 + k + 1)] = x1[2];
		X[2 * N / 4 + k + 1] = -x2[2];
		X[3 * N / 4 + k + 1] = -x3[2];
		X[0 * N / 4 + k + 1] = x0[3];
		X[1 * N / 4 + k + 1] = x1[3];
		X[N - (2 * N / 4 + k + 1)] = x2[3];
		X[N - (3 * N / 4 + k + 1)] = x3[3];
	}
	ws.destroy(std::move(Xe));
	ws.destroy(std::move(Xo));
	ws.destroy(std::move(X0));
	ws.destroy(std::move(X1));
	ws.destroy(std::move(X2));
	ws.destroy(std::move(X3));
}
