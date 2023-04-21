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

#define NTRIAL 5

void fft2simd_real(int N1, double* X, int N);

void fft_real(double* X, int N) {
	if (N <= SFFT_NMAX) {
		sfft_real(X, N);
		return;
	} else {
		static std::unordered_map<int, int> cache;
		auto iter = cache.find(N);
		if (iter == cache.end()) {
			for (int R = 4; R <= std::min(SFFT_NMAX, N); R++) {
				if (N % R == 0) {
					std::vector<double> X(N);
					for (int n = 0; n < N; n++) {
						X[n] = rand1();
					}
					fft2simd_real(R, X.data(), N);
				}
			}
			double best_time = 1e99;
			int best_R = -1;
			for (int R = 4; R <= std::min(SFFT_NMAX, N); R++) {
				if (N % R == 0) {
					std::vector<double> times;
					std::vector<double> X(N);
					for (int k = 0; k < NTRIAL; k++) {
						for (int n = 0; n < N; n++) {
							X[n] = rand1();
						}
						timer tm;
						tm.start();
						fft2simd_real(R, X.data(), N);
						tm.stop();
						times.push_back(tm.read());
					}
					std::sort(times.begin(), times.end());
					double tm = times[(times.size() + 1) / 2];
					if (tm < best_time) {
						best_time = tm;
						best_R = R;
					}
				}
			}
			cache[N] = best_R;
			iter = cache.find(N);
		}
		fft2simd_real(iter->second, X, N);
	}
}

std::vector<fft_method> possible_ffts_real(int N) {
	constexpr int FFT_NMAX = 32;
	std::vector<fft_method> ffts;
	fft_method m;
	for (m.R = 2; m.R <= std::min(SFFT_NMAX, N); m.R++) {
		if (N % m.R == 0) {
			m.type = FFT_CT;
			ffts.push_back(m);
		}
	}
	return ffts;
}

template<class T>
void fft_real(const fft_method& method, T* X, int N) {
	switch (method.type) {
	case FFT_CT:
		return fft_cooley_tukey_real(method.R, X, N);
	}
}

void fft_indices_real(const fft_method& method, int* I, int N) {
	switch (method.type) {
	case FFT_CT:
		return fft_cooley_tukey_indices_real(method.R, I, N);
	}
}

fft_method select_fft_real(int N) {
	static std::unordered_map<int, fft_method> cache;
	auto iter = cache.find(N);
	if (iter == cache.end()) {
		std::vector<fft_simd4> X(N);
		const auto tests = possible_ffts_real(N);
		const int M = tests.size();
		std::vector<double> timers(M);
		for (int m = 0; m < M; m++) {
			for (int n = 0; n < N; n++) {
				X[n] = rand1();
			}
			fft_real(tests[m], X.data(), N);
			std::vector<double> times;
			for (int n = 0; n < NTRIAL; n++) {
				timer tm;
				for (int n = 0; n < N; n++) {
					X[n] = rand1();
				}
				tm.start();
				fft_real(tests[m], X.data(), N);
				tm.stop();
				times.push_back(tm.read());
			}
			std::sort(times.begin(), times.end());
			timers[m] = times[(times.size() + 1) / 2];
		}
		fft_method best_method;
		double best_time = std::numeric_limits<double>::max();
		for (int m = 0; m < M; m++) {
			if (timers[m] < best_time) {
				best_time = timers[m];
				best_method = tests[m];
			}
		}
		//	printf("%i - %s\n", N, fft_method_string(best_method).c_str());
		cache[N] = best_method;
		iter = cache.find(N);
	}
	return iter->second;
}

void fft_real(fft_simd4* X, int N) {
	if (N <= SFFT_NMAX) {
		sfft_real(X, N);
		return;
	}
	fft_real(select_fft_real(N), X, N);
}

void fft_indices_real(int* I, int N) {
	if (N <= SFFT_NMAX) {
		return;
	}
	fft_indices_real(select_fft_real(N), I, N);
}

const std::vector<int>& fft_indices_real(int N) {
	static std::unordered_map<int, std::shared_ptr<std::vector<int>>>cache;
	auto iter = cache.find(N);
	if (iter == cache.end()) {
		std::vector<int> I(N);
		std::vector<int> J(N);
		std::iota(I.begin(), I.end(), 0);
		fft_indices_real(I.data(), N);
		cache[N] = std::make_shared<std::vector<int>>(std::move(I));
		iter = cache.find(N);
	}
	return *iter->second;
}

const std::vector<int>& fft_inv_indices_real(int N) {
	static std::unordered_map<int, std::shared_ptr<std::vector<int>>>cache;
	auto iter = cache.find(N);
	if (iter == cache.end()) {
		std::vector<int> I(N);
		std::vector<int> J(N);
		std::iota(I.begin(), I.end(), 0);
		fft_indices_real(I.data(), N);
		for (int n = 0; n < N; n++) {
			J[I[n]] = n;
		}
		cache[N] = std::make_shared<std::vector<int>>(std::move(J));
		iter = cache.find(N);
	}
	return *iter->second;
}

template<class T>
void fft_scramble1_real(T* X, int N) {
	if (N <= SFFT_NMAX) {
		return;
	}
	static std::unordered_map<int, std::shared_ptr<std::vector<int>>>cache;
	auto iter = cache.find(N);
	if (iter == cache.end()) {
		std::vector<int> I = fft_inv_indices_real(N);
		cache[N] = std::make_shared<std::vector<int>>(std::move(I));
		iter = cache.find(N);
	}
	const auto& I = *iter->second;
	fft_permute(I, X);
}

template<class T>
void fft_scramble2_real(T* X, int N) {
	if (N <= SFFT_NMAX) {
		return;
	}
	static std::unordered_map<int, std::shared_ptr<std::vector<int>>>cache;
	auto iter = cache.find(N);
	if (iter == cache.end()) {
		std::vector<int> I = fft_inv_indices_real(N);
		for (int n = 1; n < N - n; n++) {
			std::swap(I[n], I[N - n]);
		}
		cache[N] = std::make_shared<std::vector<int>>(std::move(I));
		iter = cache.find(N);
	}
	const auto& I = *iter->second;
	fft_permute(I, X);
}

void fft_scramble_real(double* X, int N) {
	fft_scramble1_real(X, N);
}

void fft_scramble_real(fft_simd4* X, int N) {
	fft_scramble1_real(X, N);
}

void fft_scramble_inv_real(double* X, int N) {
	fft_scramble2_real(X, N);
}

void fft_scramble_inv_real(fft_simd4* X, int N) {
	fft_scramble2_real(X, N);
}

void fft2simd_real(int N1, double* X, int N) {
	if (N <= SFFT_NMAX) {
		sfft_real(X, N);
		return;
	}

	const int N2 = N / N1;
	const int No2 = N / 2;
	const int N1p1o2 = (N1 + 1) / 2;
	const int N2p1o2 = (N2 + 1) / 2;
	const int N1o2 = N1 / 2;
	const int N1v = round_down(N1, SIMD_SIZE);
	const int N1voS = round_down(N1, SIMD_SIZE) / SIMD_SIZE;
	const int N1s = N1 - N1v;
	const int N2p1o2v = round_down((N2p1o2 - 1), SIMD_SIZE) + 1;

	static std::stack<std::vector<fft_simd4>> vstack;
	static std::stack<std::vector<double>> sstack;
	std::vector<fft_simd4> Yv;
	std::vector<double> Ys;
	if (vstack.size()) {
		Yv = std::move(vstack.top());
		vstack.pop();
	}
	if (sstack.size()) {
		Ys = std::move(sstack.top());
		sstack.pop();
	}
	Yv.resize(N2 * N1voS);
	Ys.resize(N2 * N1s);

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
	for (int n1 = N1v; n1 < N1; n1++) {
		for (int n2 = 0; n2 < N2; n2++) {
			Ys[(n1 - N1v) * N2 + n2] = X[N1 * n2 + n1];
		}
	}

	for (int n1 = 0; n1 < N1v; n1 += SIMD_SIZE) {
		const int o = (n1 / SIMD_SIZE) * N2;
		fft_real(Yv.data() + o, N2);
	}
	for (int n1 = N1v; n1 < N1; n1++) {
		const int o = (n1 - N1v) * N2;
		fft_real(Ys.data() + o, N2);
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
			const int ir = N2 * (n1 - N1v) + k2;
			const int ii = N2 * (n1 - N1v) - k2 + N2;
			z.real() = Ys[ir];
			z.imag() = Ys[ii];
			z *= Ws[n1 * k2];
			Ys[ir] = z.real();
			Ys[ii] = z.imag();
		}
	}
	{
		double q[N1];
		for (int n1r = 0; n1r < N1voS; n1r++) {
			for (int n1c = 0; n1c < SIMD_SIZE; n1c++) {
				const int n1 = n1r * SIMD_SIZE + n1c;
				q[n1] = Yv[N2 * n1r][n1c];
			}
		}
		for (int n1 = N1v; n1 < N1; n1++) {
			q[n1] = Ys[N2 * (n1 - N1v)];
		}
		sfft_real(q, N1);
		X[0] = q[0];
		for (int k1 = 1; k1 < N1p1o2; k1++) {
			X[0 + N2 * k1] = q[k1];
			X[N - N2 * k1] = q[N1 - k1];
		}
		if (N1 % 2 == 0) {
			X[No2] = q[N1o2];
		}
	}
	for (int k2 = 1; k2 < N2p1o2v; k2 += SIMD_SIZE) {
		complex<fft_simd4> p[N1];
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
				p[n1].real()[i] = Ys[(n1 - N1v) * N2 + k2 + i];
				p[n1].imag()[i] = Ys[(n1 - N1v) * N2 - k2 - i + N2];
			}
		}
		sfft_complex((fft_simd4*) p, N1);
		for (int k1 = 0; k1 < N1; k1++) {
			for (int i = 0; i < SIMD_SIZE; i++) {
				const int k = N2 * k1 + k2 + i;
				if (k < N - k) {
					X[k] = p[k1].real()[i];
					X[N - k] = p[k1].imag()[i];
				} else {
					X[N - k] = p[k1].real()[i];
					X[k] = -p[k1].imag()[i];
				}
			}
		}
	}
	for (int k2 = N2p1o2v; k2 < N2p1o2; k2++) {
		complex<double> p[N1];
		for (int n1r = 0; n1r < N1voS; n1r++) {
			for (int n1c = 0; n1c < SIMD_SIZE; n1c++) {
				const int n1 = n1r * SIMD_SIZE + n1c;
				p[n1].real() = Yv[n1r * N2 + k2][n1c];
				p[n1].imag() = Yv[n1r * N2 - k2 + N2][n1c];
			}
		}
		for (int n1 = N1v; n1 < N1; n1++) {
			p[n1].real() = Ys[(n1 - N1v) * N2 + k2];
			p[n1].imag() = Ys[(n1 - N1v) * N2 - k2 + N2];
		}
		sfft_complex((double*) p, N1);
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
		double q[N1];
		for (int n1r = 0; n1r < N1voS; n1r++) {
			for (int n1c = 0; n1c < SIMD_SIZE; n1c++) {
				q[n1r * SIMD_SIZE + n1c] = Yv[N2 * n1r + N2 / 2][n1c];
			}
		}
		for (int n1 = N1v; n1 < N1; n1++) {
			q[n1] = Ys[N2 * (n1 - N1v) + N2 / 2];
		}
		sfft_skew(q, N1);
		for (int k1 = 0; k1 < N1p1o2; k1++) {
			X[0 + (N2 * k1 + N2 / 2)] = q[k1];
			X[N - (N2 * k1 + N2 / 2)] = q[N1 - k1 - 1];
		}
		if (N1 % 2 == 1) {
			X[No2] = q[N1o2];
		}
	}
	vstack.push(std::move(Yv));
	sstack.push(std::move(Ys));
}

