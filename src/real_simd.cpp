#include "fftu.hpp"
#include "util.hpp"

#include <vector>
#include <unordered_map>
#include <numeric>
#include <limits>
#include <cassert>
#include <algorithm>
#include <stack>
#include <memory>
#include <cstring>

void fft2simd_real_simd_width(int N1, double* X, int N) {
	if (N <= SFFT_NMAX) {
		sfft_real(X, N);
		return;
	}
	const int N2 = N / N1;
	const int No2 = N / 2;
	const int N2o2 = N2 / 2;
	const int N1p1o2 = (N1 + 1) / 2;
	const int N2p1o2 = (N2 + 1) / 2;
	const int N1o2 = N1 / 2;
	const int N1oS = N1 / SIMD_SIZE;

	workspace<fft_simd4> vws;
	workspace<double> sws;
	workspace<complex<fft_simd4>> cvws;
	workspace<complex<double>> csws;
	auto Yv = vws.create(N);
	auto q = sws.create(N1);
	auto p = cvws.create(N1);
	auto w = csws.create(N1);
	const auto& Wv = vector_twiddles(N1, N2);

	const auto& I = fft_indices_real(N2);
	for (int n2 = 0; n2 < N2; n2++) {
		for (int n1r = 0; n1r < N1oS; n1r++) {
			for (int n1c = 0; n1c < SIMD_SIZE; n1c++) {
				const int n1 = n1r * SIMD_SIZE + n1c;
				Yv[n1r * N2 + n2][n1c] = X[I[n2] * N1 + n1];
			}
		}
	}
	for (int n1 = 0; n1 < N1; n1 += SIMD_SIZE) {
		const int o = (n1 / SIMD_SIZE) * N2;
		fft_real(Yv.data() + o, N2);
	}
	for (int k2 = 1; k2 < N2p1o2; k2++) {
		for (int n1r = 0; n1r < N1oS; n1r++) {
			complex<fft_simd4> z;
			const int ir = N2 * n1r + k2;
			const int ii = N2 * n1r - k2 + N2;
			z.real() = Yv[ir];
			z.imag() = Yv[ii];
			z *= Wv[n1r][k2];
			Yv[ir] = z.real();
			Yv[ii] = z.imag();
		}
	}
	for (int n1r = 0; n1r < N1oS; n1r++) {
		for (int n1c = 0; n1c < SIMD_SIZE; n1c++) {
			const int n1 = n1r * SIMD_SIZE + n1c;
			q[n1] = Yv[N2 * n1r][n1c];
		}
	}
	if (N1 <= SFFT_NMAX) {
		sfft_real(q.data(), N1);
	} else {
		fft_real(q.data(), N1);
	}
	X[0] = q[0];
	for (int k1 = 1; k1 < N1p1o2; k1++) {
		X[0 + N2 * k1] = q[k1];
		X[N - N2 * k1] = q[N1 - k1];
	}
	X[No2] = q[N1o2];
	for (int k2 = 1; k2 < N2p1o2; k2 += SIMD_SIZE) {
		for (int n1r = 0; n1r < N1oS; n1r++) {
			for (int n1c = 0; n1c < SIMD_SIZE; n1c++) {
				const int n1 = n1r * SIMD_SIZE + n1c;
				for (int i = 0; i < SIMD_SIZE; i++) {
					p[n1].real()[i] = Yv[n1r * N2 + k2 + i][n1c];
					p[n1].imag()[i] = Yv[n1r * N2 - k2 - i + N2][n1c];
				}
			}
		}
		if (N1 <= SFFT_NMAX) {
			sfft_complex((fft_simd4*) p.data(), N1);
		} else {
			fft_scramble(p.data(), N1);
			fft(p.data(), N1);
		}
		for (int k1 = 0; k1 < N1p1o2; k1++) {
			for (int i = 0; i < SIMD_SIZE; i++) {
				const int k = N2 * k1 + k2 + i;
				X[k] = p[k1].real()[i];
				X[N - k] = p[k1].imag()[i];
			}
		}
		for (int k1 = N1p1o2; k1 < N1; k1++) {
			for (int i = 0; i < SIMD_SIZE; i++) {
				const int k = N2 * k1 + k2 + i;
				X[N - k] = p[k1].real()[i];
				X[k] = -p[k1].imag()[i];
			}
		}
	}
	{
		for (int n1r = 0; n1r < N1oS; n1r++) {
			for (int n1c = 0; n1c < SIMD_SIZE; n1c++) {
				q[n1r * SIMD_SIZE + n1c] = Yv[N2 * n1r + N2o2][n1c];
			}
		}
		if (N1 <= SFFT_NMAX) {
			sfft_skew(q.data(), N1);
		} else {
			assert(false);
			abort();
		}
		for (int k1 = 0; k1 < N1p1o2; k1++) {
			X[N2 * k1 + N2o2] = q[k1];
			X[N - N2 * k1 - N2o2] = q[N1 - k1 - 1];
		}
		if (N1 % 2 == 1) {
			X[No2] = q[N1o2];
		}
	}
	vws.destroy(std::move(Yv));
	csws.destroy(std::move(w));
	cvws.destroy(std::move(p));
	sws.destroy(std::move(q));
}

void fft2simd_real(int N1, double* X, int N) {
	if (N <= SFFT_NMAX) {
		sfft_real(X, N);
		return;
	} else if (N1 % SIMD_SIZE == 0 && (N / N1) % SIMD_SIZE == 0) {
		fft2simd_real_simd_width(N1, X, N);
		return;
	}
	const int N2 = N / N1;
	const int No2 = N / 2;
	const int N2o2 = N2 / 2;
	const int N1p1o2 = (N1 + 1) / 2;
	const int N2p1o2 = (N2 + 1) / 2;
	const int N1o2 = N1 / 2;
	const int N1v = round_down(N1, SIMD_SIZE);
	const int N1voS = round_down(N1, SIMD_SIZE) / SIMD_SIZE;
	const int N1s = N1 - N1v;
	const int N2p1o2v = round_down((N2p1o2 - 1), SIMD_SIZE) + 1;

	workspace<fft_simd4> vws;
	workspace<double> sws;
	workspace<std::vector<double>> vsws;
	workspace<complex<fft_simd4>> cvws;
	workspace<complex<double>> csws;
	auto Yv = vws.create(N2 * N1voS);
	auto Ys = vsws.create(N1s);
	auto q = sws.create(N1);
	auto p = cvws.create(N1);
	auto w = csws.create(N1);
	for (int n = 0; n < N1s; n++) {
		Ys[n].resize(N2);
	}

	const auto& Ws = twiddles(N);
	const auto& Wv = vector_twiddles(N1, N2);

	const auto& I = fft_indices_real(N2);
	for (int n2 = 0; n2 < N2; n2++) {
		for (int n1r = 0; n1r < N1voS; n1r++) {
			for (int n1c = 0; n1c < SIMD_SIZE; n1c++) {
				const int n1 = n1r * SIMD_SIZE + n1c;
				Yv[n1r * N2 + n2][n1c] = X[I[n2] * N1 + n1];
			}
		}
	}
	for (int n2 = 0; n2 < N2; n2++) {
		for (int n1 = N1v; n1 < N1; n1++) {
			Ys[n1 - N1v][n2] = X[N1 * n2 + n1];
		}
	}
	for (int n1 = 0; n1 < N1v; n1 += SIMD_SIZE) {
		const int o = (n1 / SIMD_SIZE) * N2;
		fft_real(Yv.data() + o, N2);
	}
	for (int n1 = N1v; n1 < N1; n1++) {
		fft_real(Ys[n1 - N1v].data(), N2);
	}
	for (int k2 = 1; k2 < N2p1o2; k2++) {
		for (int n1r = 0; n1r < N1voS; n1r++) {
			complex<fft_simd4> z;
			const int ir = N2 * n1r + k2;
			const int ii = N2 * n1r - k2 + N2;
			z.real() = Yv[ir];
			z.imag() = Yv[ii];
			z *= Wv[n1r][k2];
			Yv[ir] = z.real();
			Yv[ii] = z.imag();
		}
		for (int n1 = N1v; n1 < N1; n1++) {
			complex<double> z;
			z.real() = Ys[n1 - N1v][k2];
			z.imag() = Ys[n1 - N1v][N2 - k2];
			z *= Ws[n1 * k2];
			Ys[n1 - N1v][k2] = z.real();
			Ys[n1 - N1v][N2 - k2] = z.imag();
		}
	}
	for (int n1r = 0; n1r < N1voS; n1r++) {
		for (int n1c = 0; n1c < SIMD_SIZE; n1c++) {
			const int n1 = n1r * SIMD_SIZE + n1c;
			q[n1] = Yv[N2 * n1r][n1c];
		}
	}
	for (int n1 = N1v; n1 < N1; n1++) {
		q[n1] = Ys[n1 - N1v][0];
	}
	if (N1 <= SFFT_NMAX) {
		sfft_real(q.data(), N1);
	} else {
		fft_real(q.data(), N1);
	}
	X[0] = q[0];
	for (int k1 = 1; k1 < N1p1o2; k1++) {
		X[0 + N2 * k1] = q[k1];
		X[N - N2 * k1] = q[N1 - k1];
	}
	if (N1 % 2 == 0) {
		X[No2] = q[N1o2];
	}
	for (int k2 = 1; k2 < N2p1o2v; k2 += SIMD_SIZE) {
		for (int n1r = 0; n1r < N1voS; n1r++) {
			for (int n1c = 0; n1c < SIMD_SIZE; n1c++) {
				const int n1 = n1r * SIMD_SIZE + n1c;
				for (int i = 0; i < SIMD_SIZE; i++) {
					p[n1].real()[i] = Yv[n1r * N2 + k2 + i][n1c];
					p[n1].imag()[i] = Yv[n1r * N2 - k2 - i + N2][n1c];
				}
			}
		}
		for (int n1 = N1v; n1 < N1; n1++) {
			for (int i = 0; i < SIMD_SIZE; i++) {
				p[n1].real()[i] = Ys[n1 - N1v][k2 + i];
				p[n1].imag()[i] = Ys[n1 - N1v][-k2 - i + N2];
			}
		}
		if (N1 <= SFFT_NMAX) {
			sfft_complex((fft_simd4*) p.data(), N1);
		} else {
			fft_scramble(p.data(), N1);
			fft(p.data(), N1);
		}
		for (int k1 = 0; k1 < N1p1o2; k1++) {
			for (int i = 0; i < SIMD_SIZE; i++) {
				const int k = N2 * k1 + k2 + i;
				X[k] = p[k1].real()[i];
				X[N - k] = p[k1].imag()[i];
			}
		}
		for (int k1 = N1p1o2; k1 < N1; k1++) {
			for (int i = 0; i < SIMD_SIZE; i++) {
				const int k = N2 * k1 + k2 + i;
				X[N - k] = p[k1].real()[i];
				X[k] = -p[k1].imag()[i];
			}
		}
	}
	{
		for (int k2 = N2p1o2v; k2 < N2p1o2; k2++) {
			for (int n1r = 0; n1r < N1voS; n1r++) {
				for (int n1c = 0; n1c < SIMD_SIZE; n1c++) {
					const int n1 = n1r * SIMD_SIZE + n1c;
					w[n1].real() = Yv[n1r * N2 + k2][n1c];
					w[n1].imag() = Yv[n1r * N2 - k2 + N2][n1c];
				}
			}
			for (int n1 = N1v; n1 < N1; n1++) {
				w[n1].real() = Ys[n1 - N1v][+k2];
				w[n1].imag() = Ys[n1 - N1v][-k2 + N2];
			}
			if (N1 <= SFFT_NMAX) {
				sfft_complex((double*) w.data(), N1);
			} else {
				fft(w.data(), N1);
			}
			for (int k1 = 0; k1 < N1p1o2; k1++) {
				const int k = N2 * k1 + k2;
				X[k] = w[k1].real();
				X[N - k] = w[k1].imag();
			}
			for (int k1 = N1p1o2; k1 < N1; k1++) {
				const int k = N2 * k1 + k2;
				X[N - k] = w[k1].real();
				X[k] = -w[k1].imag();
			}
		}
	}
	if (N2 % 2 == 0) {
		for (int n1r = 0; n1r < N1voS; n1r++) {
			for (int n1c = 0; n1c < SIMD_SIZE; n1c++) {
				q[n1r * SIMD_SIZE + n1c] = Yv[N2 * n1r + N2o2][n1c];
			}
		}
		for (int n1 = N1v; n1 < N1; n1++) {
			q[n1] = Ys[n1 - N1v][N2o2];
		}
		if (N1 <= SFFT_NMAX) {
			sfft_skew(q.data(), N1);
		} else {
			assert(false);
			abort();
		}
		for (int k1 = 0; k1 < N1p1o2; k1++) {
			X[N2 * k1 + N2o2] = q[k1];
			X[N - N2 * k1 - N2o2] = q[N1 - k1 - 1];
		}
		if (N1 % 2 == 1) {
			X[No2] = q[N1o2];
		}
	}
	vws.destroy(std::move(Yv));
	csws.destroy(std::move(w));
	vsws.destroy(std::move(Ys));
	cvws.destroy(std::move(p));
	sws.destroy(std::move(q));
}

