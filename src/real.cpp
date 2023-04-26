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

std::vector<fft_method> possible_start_ffts_real(int N) {
	std::vector<fft_method> ffts;
	fft_method m;
	auto fac = prime_factorization(N);
	bool found_ct = false;
	for (int R = SIMD_SIZE; R <= SFFT_NMAX; R++) {
		if (N % R == 0) {
			int N2 = N / R;
			if( N2 % 2 == 0 && R > SFFT_NMAX) {
				continue;
			}
			m.type = FFT_CT;
			m.R = R;
			ffts.push_back(m);
			found_ct = true;
		}
	}
	if (!found_ct) {
		int n = SFFT_NMAX + 1;
		while (N % n != 0 && n < N) {
			n++;
		}
		if (n < N) {
			int N2 = N / n;
			if( !(N2 % 2 == 0 && n > SFFT_NMAX)) {
				m.type = FFT_CT;
				m.R = n;
				ffts.push_back(m);
			}
		}
	}
	const std::vector<std::pair<int, int>> vfac(fac.begin(), fac.end());
	if (fac.size() > 2) {
		for (unsigned l = 1; l < fac.size(); l++) {
			auto nck = nchoosek(fac.size(), l);
			for (unsigned p = 0; p < nck.size(); p++) {
				int N1 = 1;
				for (unsigned m = 0; m < nck[p].size(); m++) {
					N1 *= std::pow(vfac[nck[p][m]].first, vfac[nck[p][m]].second);
				}
				m.R = N1;
				m.type = FFT_PFAC;
				ffts.push_back(m);
			}
		}
	}
	if (fac.size() == 1 && fac.begin()->first > SFFT_NMAX) {
		m.type = FFT_RADERS;
		ffts.push_back(m);
		m.type = FFT_RADERS_PADDED;
		ffts.push_back(m);
	}
	if (N % 2 == 0) {
	//	m.type = FFT_241;
	//	ffts.push_back(m);
	}
	return ffts;
}

void fft_real(const fft_method& m, double* X, int N) {
	if (N <= SFFT_NMAX) {
		sfft_real(X, N);
		return;
	}
	switch (m.type) {
	case FFT_CT:
		fft2simd_real(m.R, X, N);
		break;
	case FFT_RADERS:
		fft_raders_real(X, N, false);
		break;
	case FFT_RADERS_PADDED:
		fft_raders_real(X, N, true);
		break;
	case FFT_241:
		fft_twoforone_real(X, N);
		break;
	case FFT_PFAC:
		fft_raders_prime_factor_real(m.R, X, N);
		break;
	default:
		assert(false);
		abort();
	}
}

fft_method select_start_fft_real(int N) {
	static std::unordered_map<int, fft_method> cache;
	auto iter = cache.find(N);
	if (iter == cache.end()) {
		std::vector<double> X(N);
		const auto tests = possible_start_ffts_real(N);
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
		cache[N] = best_method;
		iter = cache.find(N);
		//printf("%i %s\n", N, fft_method_string(best_method).c_str());
	}
	return iter->second;
}

void fft_real(double* X, int N) {
	if (N <= SFFT_NMAX) {
		sfft_real(X, N);
		return;
	}
	const auto meth = select_start_fft_real(N);
	fft_real(meth, X, N);
}

std::vector<fft_method> possible_ffts_real(int N) {
	std::vector<fft_method> ffts;
	fft_method m;
	if (N % 4 == 0) {
		for (m.R = 4; m.R <= std::min(SFFT_NMAX, N); m.R += 4) {
			if (N % m.R == 0) {
				m.type = FFT_SPLIT;
				ffts.push_back(m);
			}
		}
	}
	auto pfac = prime_factorization(N);
	if (pfac.begin()->first <= SFFT_NMAX) {
		for (m.R = 2; m.R <= std::min(SFFT_NMAX, N); m.R++) {
			if (N % m.R == 0) {
			//	m.type = FFT_CONJ;
			//	ffts.push_back(m);
				m.type = FFT_CT;
				ffts.push_back(m);
			}
		}
	}
	if (pfac.size() == 1 && pfac.begin()->second == 1 && pfac.begin()->first > SFFT_NMAX) {
		m.type = FFT_RADERS;
		ffts.push_back(m);
	}
	if (ffts.size() == 0) {
		m.R = SFFT_NMAX + 1;
		while (N % m.R != 0) {
			m.R++;
		}
		m.type = FFT_CT;
		ffts.push_back(m);
	}
	assert(ffts.size());
	return ffts;
}

void fft_real(const fft_method& method, fft_simd4* X, int N) {
	switch (method.type) {
	case FFT_CT:
		fft_cooley_tukey_real(method.R, X, N);
		break;
	case FFT_CONJ:
		fft_conjugate_real(method.R, X, N);
		break;
	case FFT_RADERS:
		fft_raders_real(X, N);
		break;
	case FFT_SPLIT:
		fft_split_real(method.R, X, N);
		break;
	case FFT_241:
		fft_twoforone_real(X, N);
		break;
	default:
		assert(false);
		abort();
	}
}

void fft_raders_indices_real(int* I, int N);

void fft_indices_real(const fft_method& method, int* I, int N) {
	switch (method.type) {
	case FFT_CT:
		return fft_cooley_tukey_indices_real(method.R, I, N);
	case FFT_CONJ:
		return fft_conjugate_indices_real(method.R, I, N);
	case FFT_RADERS:
		return fft_raders_indices_real(I, N);
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
	//	printf( "%i : %s\n", N, fft_method_string(best_method).c_str());
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

