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
			for (int R = 4; R <= std::min(N, SFFT_NMAX); R++) {
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
			for (int R = 4; R <= std::min(N, SFFT_NMAX); R++) {
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
		if( iter->second < 0 ) {
			fft_raders_real(X, N, false);
		} else {
			fft2simd_real(iter->second, X, N);
		}
	}
}

std::vector<fft_method> possible_ffts_real(int N) {
	std::vector<fft_method> ffts;
	fft_method m;
	if (N % 4 == 0) {
		m.R = 4;
		m.type = FFT_SPLIT;
		ffts.push_back(m);
		if (N % 2 == 0) {
			m.type = FFT_241;
			ffts.push_back(m);
		}
	}
	auto pfac = prime_factorization(N);
	if (pfac.rbegin()->first <= SFFT_NMAX) {
		for (m.R = 2; m.R <= std::min(SFFT_NMAX, N); m.R++) {
			if (N % m.R == 0) {
				m.type = FFT_CT;
				ffts.push_back(m);
			}
		}
	}
	return ffts;
}

template<class T>
void fft_real(const fft_method& method, T* X, int N) {
	switch (method.type) {
	case FFT_CT:
		fft_cooley_tukey_real(method.R, X, N);
		break;
	case FFT_SPLIT:
		fft_split_real(method.R, X, N);
		break;
	case FFT_241:
		fft_twoforone_real(X, N);
		break;
	}
}

void fft_indices_real(const fft_method& method, int* I, int N) {
	switch (method.type) {
	case FFT_CT:
		return fft_cooley_tukey_indices_real(method.R, I, N);
	case FFT_SPLIT:
		return fft_split_real_indices(method.R, I, N);
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
//		printf("%i - %s\n", N, fft_method_string(best_method).c_str());
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

