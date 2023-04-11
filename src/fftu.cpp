#include "fftu.hpp"
#include "util.hpp"
#include "fft.hpp"

#include <vector>
#include <unordered_map>
#include <numeric>
#include <limits>
#include <cassert>
#include <cstring>

#define FFT_SPLIT 3
#define FFT_CT 1
#define FFT_6STEP 2

struct fft_method {
	int type;
	int R;
};

void fft_scramble(fft_method method, complex<fft_simd4>* X, int N);

std::string fft_method_string(const fft_method method) {
	std::string str;
	switch (method.type) {
	case FFT_SPLIT:
		str = "split-" + std::to_string(method.R);
		break;
	case FFT_CT:
		str = "cooley-tukey-" + std::to_string(method.R);
		break;
	case FFT_6STEP:
		str = "six-step";
		break;
	default:
		assert(false);
	}
	return str;
}

void fft(complex<double>* X, int N) {
	static std::vector<complex<fft_simd4>> Y;
	Y.resize(N / SIMD_SIZE);
	if (N <= FFT_NMAX) {
		switch (N) {
		case 2:
			return fft_complex_exe_2((double*) X);
		case 4:
			return fft_complex_exe_4((double*) X);
		case 8:
			return fft_complex_exe_8((double*) X);
		case 16:
			return fft_complex_exe_16((double*) X);
		case 32:
			return fft_complex_exe_32((double*) X);
		}
		return;
	}
	constexpr int N1 = 16;
	const int N2 = N / N1;
	const auto& I = fft_indices(N2);
	std::array<complex<fft_simd4>, N1> z;
	const auto& W = vector_twiddles(N1, N2);
	for (int n2 = 0; n2 < N2; n2++) {
		for (int n1 = 0; n1 < N1; n1++) {
			Y[(n1 / SIMD_SIZE) * N2 + I[n2]].real()[n1 % SIMD_SIZE] = X[N1 * n2 + n1].real();
			Y[(n1 / SIMD_SIZE) * N2 + I[n2]].imag()[n1 % SIMD_SIZE] = X[N1 * n2 + n1].imag();
		}
	}
	for (int n1 = 0; n1 < N1; n1 += SIMD_SIZE) {
		const int o = (n1 / SIMD_SIZE) * N2;
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
		fft_complex_exe_16((fft_simd4*) z.data());
		for (int k1 = 0; k1 < N1; k1++) {
			for (int i = 0; i < SIMD_SIZE; i++) {
				X[k1 * N2 + k2 + i].real() = z[k1].real()[i];
				X[k1 * N2 + k2 + i].imag() = z[k1].imag()[i];
			}
		}
	}
	Y.resize(0);
}

std::vector<fft_method> possible_ffts(int N) {
	std::vector<fft_method> ffts;
	fft_method m;
	const int logN = ilogb(N);
	if (logN % 2 == 0) {
		m.type = FFT_6STEP;
		ffts.push_back(m);
	}
	m.type = FFT_SPLIT;
	for (m.R = 4; m.R <= FFT_NMAX; m.R *= 2) {
		ffts.push_back(m);
	}
	m.type = FFT_CT;
	for (m.R = 2; m.R <= FFT_NMAX; m.R *= 2) {
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
	case FFT_6STEP:
		return fft_six_step(X, N);
	}
}

void fft_indices(const fft_method& method, int* I, int N) {
	switch (method.type) {
	case FFT_SPLIT:
		return fft_split_indices(method.R, I, N);
	case FFT_CT:
		return fft_cooley_tukey_indices(method.R, I, N);
	case FFT_6STEP:
		return;
	}
}

fft_method fft_select(int N) {
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
			int ntests = 20;
			timer tm;
			double best = 1e99;
			double worst = -1e99;
			for (int l = 0; l < ntests; l++) {
				tm.start();
				fft(tests[m], X.data(), N);
				tm.stop();
				if (tm.read() < best) {
					best = tm.read();
				}
				if (tm.read() > worst) {
					worst = tm.read();
				}
				timers[m] += tm.read();
			}
			timers[m] -= best + worst;
		}
		fft_method best_method;
		double best_time = std::numeric_limits<double>::max();
		for (int m = 0; m < M; m++) {
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
	if (N <= 16) {
		switch (N) {
		case 2:
			return fft_complex_exe_2((fft_simd4*) X);
		case 4:
			return fft_complex_exe_4((fft_simd4*) X);
		case 8:
			return fft_complex_exe_8((fft_simd4*) X);
		case 16:
			return fft_complex_exe_16((fft_simd4*) X);
		case 32:
			return fft_complex_exe_32((fft_simd4*) X);
		}
		return;
	}
	fft(fft_select(N), X, N);
}

void fft_indices(int* I, int N) {
	if (N <= FFT_NMAX) {
		return;
	}
	fft_indices(fft_select(N), I, N);
}

const std::vector<int>& fft_indices(int N) {
	static std::unordered_map<int, std::vector<int>> cache;
	auto iter = cache.find(N);
	if (iter == cache.end()) {
		std::vector<int> I(N);
		std::vector<int> J(N);
		std::iota(I.begin(), I.end(), 0);
		fft_indices(fft_select(N), I.data(), N);
		for (int n = 0; n < N; n++) {
			J[I[n]] = n;
		}
		cache[N] = std::move(J);
		iter = cache.find(N);
	}
	return iter->second;
}

void fft_scramble(complex<fft_simd4>* X, int N) {
	fft_scramble(fft_select(N), X, N);
}

void fft_scramble(fft_method method, complex<fft_simd4>* X, int N) {
	if (N <= FFT_NMAX) {
		return;
	}
	static std::unordered_map<int, std::vector<int>> cache;
	auto iter = cache.find(N);
	if (iter == cache.end()) {
		std::vector<int> I(N);
		std::vector<int> J(N);
		std::iota(I.begin(), I.end(), 0);
		fft_indices(method, I.data(), N);
		for (int n = 0; n < N; n++) {
			J[I[n]] = n;
		}
		cache[N] = std::move(J);
		iter = cache.find(N);
	}
	const auto& I = iter->second;
	/*auto Y = create_scratch(N);
	 for( int n = 0; n < N; n++) {
	 Y[I[n]] = X[n];
	 }
	 std::memcpy(X, Y.data(), sizeof(complex<fft_simd4>) * N);*/
	static std::vector<bool> flag;
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
	flag.resize(0);
//	destroy_scratch(std::move(Y));
}

