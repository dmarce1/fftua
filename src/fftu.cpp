#include "fftu.hpp"
#include "util.hpp"
#include "fft.hpp"

#include <vector>
#include <unordered_map>
#include <numeric>
#include <limits>
#include <cassert>

#define FFT_SPLIT 0
#define FFT_CT 1

struct fft_method {
	int type;
	int R;
};

std::string fft_method_string(const fft_method method) {
	std::string str;
	switch (method.type) {
	case FFT_SPLIT:
		str = "split-" + std::to_string(method.R);
		break;
	case FFT_CT:
		str = "cooley-tukey-" + std::to_string(method.R);
		break;
	default:
		assert(false);
	}
	return str;
}

void fft(complex<double>* X, int N) {
	static std::vector<complex<fft_simd4>> Y;
	Y.resize(N / SIMD_SIZE);
	if (N < FFT_NMAX) {
		fft_complex((double*) X, N);
		return;
	}
	constexpr int N1 = 16;
	const int N2 = N / N1;
	std::array<complex<fft_simd4>, N1> z;
	const auto& W = vector_twiddles(N1, N2);
	for (int n2 = 0; n2 < N2; n2++) {
		for (int n1 = 0; n1 < N1; n1++) {
			Y[(n1 / SIMD_SIZE) * N2 + n2].real()[n1 % SIMD_SIZE] = X[N1 * n2 + n1].real();
			Y[(n1 / SIMD_SIZE) * N2 + n2].imag()[n1 % SIMD_SIZE] = X[N1 * n2 + n1].imag();
		}
	}
	for (int n1 = 0; n1 < N1; n1 += SIMD_SIZE) {
		const int o = (n1 / SIMD_SIZE) * N2;
		fft_scramble((complex<fft_simd4>*) Y.data() + o, N2);
		fft((complex<fft_simd4>*) Y.data() + o, N2);
	}
	for (int k2 = 0; k2 < N2; k2++) {
		for (int n1 = 0; n1 < N1; n1 += SIMD_SIZE) {
			Y[(n1 / SIMD_SIZE) * N2 + k2] *= W[n1 / SIMD_SIZE][k2];
		}
	}
	for (int k2 = 0; k2 < N2; k2 += SIMD_SIZE) {
		for (int n1 = 0; n1 < N1; n1++) {
			for (int i = 0; i < SIMD_SIZE; i++) {
				z[n1].real()[i] = Y[(n1 / SIMD_SIZE) * N2 + k2 + i].real()[n1 % SIMD_SIZE];
				z[n1].imag()[i] = Y[(n1 / SIMD_SIZE) * N2 + k2 + i].imag()[n1 % SIMD_SIZE];
			}
		}
		fft_complex_simd4((fft_simd4*) z.data(), N1);
		for (int k1 = 0; k1 < N1; k1++) {
			for (int i = 0; i < SIMD_SIZE; i++) {
				X[k1 * N2 + k2 + i].real() = z[k1].real()[i];
				X[k1 * N2 + k2 + i].imag() = z[k1].imag()[i];
			}
		}
	}
}

std::vector<fft_method> possible_ffts(int N) {
	std::vector<fft_method> ffts;
	fft_method m;
	m.type = FFT_SPLIT;
	for (m.R = 4; m.R <= FFT_NMAX; m.R *= 2) {
		ffts.push_back(m);
	}
	m.type = FFT_CT;
	for (m.R = 2; m.R < sqrt(N); m.R *= 2) {
		ffts.push_back(m);
	}
	return ffts;
}

void fft(const fft_method& method, complex<fft_simd4>* X, int N) {
	switch (method.type) {
	case FFT_SPLIT:
		return fft_split(method.R, X, N);
	case FFT_CT:
		return fft_cooley_tukey(method.R, X, N);
	}
}

void fft_indices(const fft_method& method, int* I, int N) {
	switch (method.type) {
	case FFT_SPLIT:
		return fft_split_indices(method.R, I, N);
	case FFT_CT:
		return fft_cooley_tukey_indices(method.R, I, N);
	}
}

fft_method select_fft(int N) {
	static std::unordered_map<int, fft_method> cache;
	auto iter = cache.find(N);
	if (iter == cache.end()) {
		std::vector<complex<fft_simd4>> X(N);
		for (int n = 0; n < N; n++) {
			X[n].real() = rand1();
			X[n].imag() = rand1();
		}
		const auto tests = possible_ffts(N);
		const int M = tests.size();
		std::vector<double> timers(M);
		for (int m = 0; m < M; m++) {
			fft(tests[m], X.data(), N);
			timer t1, t2, t3;
			t1.start();
			fft(tests[m], X.data(), N);
			t1.stop();
			t2.start();
			fft(tests[m], X.data(), N);
			t2.stop();
			t3.start();
			fft(tests[m], X.data(), N);
			t3.stop();
			if ((t3.read() - t1.read()) * (t1.read() - t2.read()) <= 0.0) {
				timers[m] = 0.1 * t1.read();
			} else if ((t3.read() - t2.read()) * (t2.read() - t1.read()) <= 0.0) {
				timers[m] = 0.1 * t2.read();
			} else if ((t1.read() - t3.read()) * (t3.read() - t2.read()) <= 0.0) {
				timers[m] = 0.1 * t3.read();
			}
			timers[m] += 0.9 * (t1.read() + t2.read() + t3.read());
		}
		fft_method best_method;
		double best_time = std::numeric_limits<double>::max();
		for (int m = 0; m < M; m++) {
			//		printf("%s - %e\n", fft_method_string(tests[m]).c_str(), timers[m]);
			if (timers[m] < best_time) {
				best_time = timers[m];
				best_method = tests[m];
			}
		}
		printf("%s\n\n", fft_method_string(best_method).c_str());
		cache[N] = best_method;
		iter = cache.find(N);
	}
	return iter->second;
}

void fft(complex<fft_simd4>* X, int N) {
	if (N <= FFT_NMAX) {
		fft_complex_simd4((fft_simd4*) X, N);
		return;
	}
	fft(select_fft(N), X, N);
}

void fft_indices(int* I, int N) {
	if (N <= FFT_NMAX) {
		return;
	}
	fft_indices(select_fft(N), I, N);
}

void fft_scramble(complex<fft_simd4>* X, int N) {
	if (N <= FFT_NMAX) {
		return;
	}
	static std::unordered_map<int, std::vector<int>> cache;
	auto iter = cache.find(N);
	if (iter == cache.end()) {
		std::vector<int> I(N);
		std::vector<int> J(N);
		std::iota(I.begin(), I.end(), 0);
		fft_indices(I.data(), N);
		for (int n = 0; n < N; n++) {
			J[I[n]] = n;
		}
		cache[N] = std::move(J);
		iter = cache.find(N);
	}
	const auto& I = iter->second;
	static std::vector<bool> flag;
	flag.resize(0);
	flag.resize(N, false);
	for (int n = 0; n < N; n++) {
		if (!flag[n]) {
			flag[n] = true;
			complex<fft_simd4> tmp = X[n];
			int m = I[n];
			while (m != n) {
				flag[m] = true;
				std::swap(tmp, X[m]);
				m = I[m];
			}
			X[n] = tmp;
		}
	}
}

