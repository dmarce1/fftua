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
#define FFT_6 2

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
	case FFT_6:
		str = "6-step";
		break;
	default:
		assert(false);
	}
	return str;
}

template<int N1>
void fft(complex<double>* X, int N) {
	static std::vector<complex<fft_simd4>> Y;
	static std::vector<complex<fft_simd4>> Z;
	if (N <= FFT_NMAX) {
		switch (N) {
		case 2:
			return fft_complex_2((double*) X);
		case 4:
			return fft_complex_4((double*) X);
		case 8:
			return fft_complex_8((double*) X);
		case 16:
			return fft_complex_16((double*) X);
		case 32:
			return fft_complex_32((double*) X);
		case 64:
			return fft_complex_64((double*) X);
		}
		return;
	}
	Y.resize(N / SIMD_SIZE);
	Z.resize(N / SIMD_SIZE);
	const int N2 = N / N1;
	std::array<complex<fft_simd4>, N1> z;
	const auto& I = fft_indices(N2);
	const auto& W = vector_twiddles(N1, N2);
	for (int n2 = 0; n2 < N2; n2++) {
		for (int n1c = 0; n1c < SIMD_SIZE; n1c++) {
			for (int n1r = 0; n1r < N1 / SIMD_SIZE; n1r++) {
				const int n1 = n1r * SIMD_SIZE + n1c;
				Y[n1r * N2 + I[n2]].real()[n1c] = X[N1 * n2 + n1].real();
				Y[n1r * N2 + I[n2]].imag()[n1c] = X[N1 * n2 + n1].imag();
			}
		}
	}
	for (int n1 = 0; n1 < N1; n1 += SIMD_SIZE) {
		const int o = (n1 / SIMD_SIZE) * N2;
		fft((complex<fft_simd4>*) Y.data() + o, (complex<fft_simd4>*) Z.data() + o, N2);
	}
	for (int k2 = 0; k2 < N2; k2++) {
		for (int n1 = 0; n1 < N1; n1 += SIMD_SIZE) {
			Y[(n1 / SIMD_SIZE) * N2 + k2] *= W[n1 / SIMD_SIZE][k2];
		}
	}
	for (int k2 = 0; k2 < N2; k2 += SIMD_SIZE) {
		for (int n1r = 0; n1r < N1 / SIMD_SIZE; n1r++) {
			for (int n1c = 0; n1c < SIMD_SIZE; n1c++) {
				const int n1 = n1r * SIMD_SIZE + n1c;
				for (int i = 0; i < SIMD_SIZE; i++) {
					z[n1].real()[i] = Y[(n1r) * N2 + k2 + i].real()[n1c];
					z[n1].imag()[i] = Y[(n1r) * N2 + k2 + i].imag()[n1c];
				}
			}
		}
		switch (N1) {
		case 4:
			fft_complex_4((fft_simd4*) z.data());
			break;
		case 8:
			fft_complex_8((fft_simd4*) z.data());
			break;
		case 16:
			fft_complex_16((fft_simd4*) z.data());
			break;
		case 32:
			fft_complex_32((fft_simd4*) z.data());
			break;
		case 64:
			fft_complex_64((fft_simd4*) z.data());
			break;
		}
		for (int k1 = 0; k1 < N1; k1++) {
			for (int i = 0; i < SIMD_SIZE; i++) {
				X[k1 * N2 + k2 + i].real() = z[k1].real()[i];
				X[k1 * N2 + k2 + i].imag() = z[k1].imag()[i];
			}
		}
	}
}

void fft(complex<double>* X, int N) {
	static std::unordered_map<int, int> cache;
	auto iter = cache.find(N);
	if (iter == cache.end()) {
		std::vector<complex<double>> Y(N);
		auto* X = Y.data();
		int best;
		fft<4>(X, N);
		fft<8>(X, N);
		fft<16>(X, N);
		if (N >= 128) {
			fft<32>(X, N);
		}
		if (N >= 256) {
			fft<64>(X, N);
		}
		timer tm;
		double best_time = 1e99;
		tm.start();
		for (int n = 0; n < 10; n++) {
			fft<4>(X, N);
		}
		tm.stop();
		best_time = tm.read();
		best = 4;
		tm.reset();
		tm.start();
		for (int n = 0; n < 10; n++) {
			fft<8>(X, N);
		}
		tm.stop();
		if (best_time > tm.read()) {
			best_time = tm.read();
			best = 8;
		}
		tm.reset();
		tm.start();
		for (int n = 0; n < 10; n++) {
			fft<16>(X, N);
		}
		tm.stop();
		if (best_time > tm.read()) {
			best_time = tm.read();
			best = 16;
		}
		tm.reset();
		if (N >= 128) {
			tm.start();
			for (int n = 0; n < 10; n++) {
				fft<32>(X, N);
			}
			tm.stop();
			if (best_time > tm.read()) {
				best_time = tm.read();
				best = 32;
			}
			tm.reset();
		}
		if (N >= 256) {
			tm.start();
			for (int n = 0; n < 10; n++) {
				fft<64>(X, N);
			}
			tm.stop();
			if (best_time > tm.read()) {
				best_time = tm.read();
				best = 64;
			}
			tm.reset();
		}
		cache[N] = best;
		iter = cache.find(N);
		printf("simd width = %i\n", cache[N]);
	}
	switch (cache[N]) {
	case 4:
		return fft<4>(X, N);
	case 8:
		return fft<8>(X, N);
	case 16:
		return fft<16>(X, N);
	case 32:
		return fft<32>(X, N);
	case 64:
		return fft<64>(X, N);
	}
}

std::vector<fft_method> possible_ffts(int N) {
	std::vector<fft_method> ffts;
	fft_method m;
	if (ilogb(N) % 2 == 0) {
		m.type = FFT_6;
		ffts.push_back(m);
	} //else {
	m.type = FFT_SPLIT;
	for (m.R = 4; m.R <= FFT_NMAX; m.R *= 2) {
		ffts.push_back(m);
	}
	m.type = FFT_CT;
	for (m.R = 2; m.R <= FFT_NMAX; m.R *= 2) {
		ffts.push_back(m);
	}
//	}
	return ffts;
}

void fft(const fft_method& method, complex<fft_simd4>* X, complex<fft_simd4>* Y, int N) {
	switch (method.type) {
	case FFT_SPLIT:
		return fft_split(method.R, X, Y, N);
	case FFT_CT:
		return fft_cooley_tukey(method.R, X, Y, N);
	case FFT_6:
		return fft_six_step(X, Y, N);
	}
}

void fft_indices(const fft_method& method, int* I, int N) {
	switch (method.type) {
	case FFT_SPLIT:
		return fft_split_indices(method.R, I, N);
	case FFT_CT:
		return fft_cooley_tukey_indices(method.R, I, N);
	case FFT_6:
		return fft_six_step_indices(I, N);
	}
}

fft_method select_fft(int N) {
	static std::unordered_map<int, fft_method> cache;
	auto iter = cache.find(N);
	if (iter == cache.end()) {
		std::vector<complex<fft_simd4>> X(N);
		std::vector<complex<fft_simd4>> Y(N);
		for (int n = 0; n < N; n++) {
			X[n].real() = rand1();
			X[n].imag() = rand1();
		}
		const auto tests = possible_ffts(N);
		const int M = tests.size();
		std::vector<double> timers(M);
		for (int m = 0; m < M; m++) {
			fft(tests[m], X.data(), Y.data(), N);
			double best = -1.0;
			double worst = 1e99;
			for (int n = 0; n < 20; n++) {
				timer tm;
				tm.start();
				fft(tests[m], X.data(), Y.data(), N);
				tm.stop();
				auto t = tm.read();
				timers[m] += t;
				tm.reset();
				best = std::max(t, best);
				worst = std::min(t, worst);
			}
			timers[m] -= best + worst;
		}
		fft_method best_method;
		double best_time = std::numeric_limits<double>::max();
		for (int m = 0; m < M; m++) {
			//	printf("%s - %e\n", fft_method_string(tests[m]).c_str(), timers[m]);
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

void fft(complex<fft_simd4>* X, complex<fft_simd4>* Y, int N) {
	if (N <= FFT_NMAX) {
		switch (N) {
		case 2:
			return fft_complex_2((fft_simd4*) X);
		case 4:
			return fft_complex_4((fft_simd4*) X);
		case 8:
			return fft_complex_8((fft_simd4*) X);
		case 16:
			return fft_complex_16((fft_simd4*) X);
		case 32:
			return fft_complex_32((fft_simd4*) X);
		case 64:
			return fft_complex_64((fft_simd4*) X);
		}
		return;
	}
	fft(select_fft(N), X, Y, N);
}

void fft_indices(int* I, int N) {
	if (N <= FFT_NMAX) {
		return;
	}
	fft_indices(select_fft(N), I, N);
}

const std::vector<int>& fft_indices(int N) {
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
	return iter->second;
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
		for (int n = 0; n < N; n++) {
			I[J[n]] = n;
		}
		cache[N] = std::move(I);
		iter = cache.find(N);
	}
	const auto& I = iter->second;
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
}

