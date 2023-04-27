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

#define NTRIAL 5

fft_method select_fft(int N);

std::string fft_method_string(const fft_method method) {
	std::string str;
	switch (method.type) {
	case FFT_SPLIT:
		str = "split-" + std::to_string(method.R);
		break;
	case FFT_SPLIT_CONJ:
		str = "split-conjugate-" + std::to_string(method.R);
		break;
	case FFT_CT:
		str = "cooley-tukey-" + std::to_string(method.R);
		break;
	case FFT_CONJ:
		str = "conjugate-" + std::to_string(method.R);
		break;
	case FFT_6:
		str = "6-step";
		break;
	case FFT_RADERS:
		str = "raders";
		break;
	case FFT_RADERS_PADDED:
		str = "raders-padded";
		break;
	case FFT_PFAC:
		str = "prime-factor-" + std::to_string(method.R);
		;
		break;
	case FFT_241:
		str = "twoforone-real";
		break;
	default:
		assert(false);
	}
	return str;
}

void fft2simd_raders(int N1, complex<double>* X, int N);

void fft2simd(int N1, complex<double>* X, int N) {
	if (N1 > SFFT_NMAX) {
		fft2simd_raders(N1, X, N);
		return;
	}
	if (N <= SFFT_NMAX) {
		sfft_complex((double*) X, N);
		return;
	}
	const int N2 = N / N1;
	const int N1v = round_down(N1, SIMD_SIZE);
	const int N1voS = round_down(N1, SIMD_SIZE) / SIMD_SIZE;
	const int N1s = N1 - N1v;
	const int N2v = round_down(N2, SIMD_SIZE);
	workspace<complex<fft_simd4>> vws;
	workspace<complex<double>> sws;
	auto Yv = vws.create(N2 * N1voS);
	auto Ys = sws.create(N2 * N1s);
	complex<fft_simd4> z[N1];
	select_fft(N2);
	const auto& I = fft_indices(N2);
	const auto& Wv = vector_twiddles(N1, N2);
	const auto& Ws = twiddles(N);
	for (int n2 = 0; n2 < N2; n2++) {
		for (int n1r = 0; n1r < N1voS; n1r++) {
			for (int n1c = 0; n1c < SIMD_SIZE; n1c++) {
				const int n1 = n1r * SIMD_SIZE + n1c;
				const int to = n1r * N2 + n2;
				const int from = N1 * I[n2] + n1;
				Yv[to].real()[n1c] = X[from].real();
				Yv[to].imag()[n1c] = X[from].imag();
			}
		}
	}
	for (int n1 = N1v; n1 < N1; n1++) {
		for (int n2 = 0; n2 < N2; n2++) {
			Ys[(n1 - N1v) * N2 + n2] = X[N1 * n2 + n1];
		}
	}
	if (N1v) {
		for (int n1 = 0; n1 < N1v; n1 += SIMD_SIZE) {
			const int o = (n1 / SIMD_SIZE) * N2;
			fft(Yv.data() + o, N2);
		}
	}
	if (N1 - N1v) {
		for (int n1 = N1v; n1 < N1; n1++) {
			const int o = (n1 - N1v) * N2;
			fft(Ys.data() + o, N2);
		}
	}
	for (int k2 = 0; k2 < N2; k2++) {
		for (int n1 = 0; n1 < N1v; n1 += SIMD_SIZE) {
			Yv[(n1 / SIMD_SIZE) * N2 + k2] *= Wv[n1 / SIMD_SIZE][k2];
		}
		for (int n1 = N1v; n1 < N1; n1++) {
			Ys[(n1 - N1v) * N2 + k2] *= Ws[n1 * k2];
		}
	}
	for (int k2 = 0; k2 < N2v; k2 += SIMD_SIZE) {
		for (int n1r = 0; n1r < N1voS; n1r++) {
			for (int n1c = 0; n1c < SIMD_SIZE; n1c++) {
				const int n1 = n1r * SIMD_SIZE + n1c;
				for (int i = 0; i < SIMD_SIZE; i++) {
					z[n1].real()[i] = Yv[n1r * N2 + k2 + i].real()[n1c];
					z[n1].imag()[i] = Yv[n1r * N2 + k2 + i].imag()[n1c];
				}
			}
		}
		for (int n1 = N1v; n1 < N1; n1++) {
			for (int i = 0; i < SIMD_SIZE; i++) {
				z[n1].real()[i] = Ys[(n1 - N1v) * N2 + k2 + i].real();
				z[n1].imag()[i] = Ys[(n1 - N1v) * N2 + k2 + i].imag();
			}
		}
		if (N1 <= SFFT_NMAX) {
			sfft_complex((fft_simd4*) z, N1);
		} else {
			fft_scramble(z, N1);
			fft(z, N1);
		}
		for (int k1 = 0; k1 < N1; k1++) {
			for (int i = 0; i < SIMD_SIZE; i++) {
				X[k1 * N2 + k2 + i].real() = z[k1].real()[i];
				X[k1 * N2 + k2 + i].imag() = z[k1].imag()[i];
			}
		}
	}
	for (int k2 = N2v; k2 < N2; k2++) {
		complex<double> z[N1];
		for (int n1r = 0; n1r < N1voS; n1r++) {
			for (int n1c = 0; n1c < SIMD_SIZE; n1c++) {
				const int n1 = n1r * SIMD_SIZE + n1c;
				z[n1].real() = Yv[n1r * N2 + k2].real()[n1c];
				z[n1].imag() = Yv[n1r * N2 + k2].imag()[n1c];
			}
		}
		for (int n1 = N1v; n1 < N1; n1++) {
			z[n1] = Ys[(n1 - N1v) * N2 + k2];
		}
		if (N1 <= SFFT_NMAX) {
			sfft_complex((double*) z, N1);
		} else {
			fft(z, N1);
		}
		for (int k1 = 0; k1 < N1; k1++) {
			X[k1 * N2 + k2] = z[k1];
		}
	}
	vws.destroy(std::move(Yv));
	sws.destroy(std::move(Ys));
}

void fft2simd_raders(int N1, complex<double>* X, int N) {
	if (N <= SFFT_NMAX) {
		sfft_complex((double*) X, N);
		return;
	}
	const int N1v = round_down(N1, SIMD_SIZE);
	const int N1voS = round_down(N1, SIMD_SIZE) / SIMD_SIZE;
	const int N1s = N1 - N1v;
	const int N2 = N / N1;
	const int N2v = round_down(N2, SIMD_SIZE);
	workspace<complex<fft_simd4>> vws;
	workspace<complex<double>> sws;
	auto Yv = vws.create(N2 * N1voS);
	auto Ys = sws.create(N2 * N1s);
	auto zv = vws.create(N1);
	auto zs = sws.create(N1);
	const auto& I = fft_indices(N2);
	const auto& Wv = vector_twiddles(N1, N2);
	const auto& Ws = twiddles(N);
	for (int n2 = 0; n2 < N2; n2++) {
		for (int n1r = 0; n1r < N1voS; n1r++) {
			for (int n1c = 0; n1c < SIMD_SIZE; n1c++) {
				const int n1 = n1r * SIMD_SIZE + n1c;
				const int to = n1r * N2 + n2;
				const int from = N1 * I[n2] + n1;
				Yv[to].real()[n1c] = X[from].real();
				Yv[to].imag()[n1c] = X[from].imag();
			}
		}
	}
	for (int n1 = N1v; n1 < N1; n1++) {
		for (int n2 = 0; n2 < N2; n2++) {
			Ys[(n1 - N1v) * N2 + n2] = X[N1 * n2 + n1];
		}
	}
	if (N1v) {
		for (int n1 = 0; n1 < N1v; n1 += SIMD_SIZE) {
			const int o = (n1 / SIMD_SIZE) * N2;
			fft(Yv.data() + o, N2);
		}
	}
	if (N1 - N1v) {
		for (int n1 = N1v; n1 < N1; n1++) {
			const int o = (n1 - N1v) * N2;
			fft(Ys.data() + o, N2);
		}
	}
	for (int k2 = 0; k2 < N2; k2++) {
		for (int n1 = 0; n1 < N1v; n1 += SIMD_SIZE) {
			Yv[(n1 / SIMD_SIZE) * N2 + k2] *= Wv[n1 / SIMD_SIZE][k2];
		}
		for (int n1 = N1v; n1 < N1; n1++) {
			Ys[(n1 - N1v) * N2 + k2] *= Ws[n1 * k2];
		}
	}
	for (int k2 = 0; k2 < N2v; k2 += SIMD_SIZE) {
		for (int n1r = 0; n1r < N1voS; n1r++) {
			for (int n1c = 0; n1c < SIMD_SIZE; n1c++) {
				const int n1 = n1r * SIMD_SIZE + n1c;
				for (int i = 0; i < SIMD_SIZE; i++) {
					zv[n1].real()[i] = Yv[n1r * N2 + k2 + i].real()[n1c];
					zv[n1].imag()[i] = Yv[n1r * N2 + k2 + i].imag()[n1c];
				}
			}
		}
		for (int n1 = N1v; n1 < N1; n1++) {
			for (int i = 0; i < SIMD_SIZE; i++) {
				zv[n1].real()[i] = Ys[(n1 - N1v) * N2 + k2 + i].real();
				zv[n1].imag()[i] = Ys[(n1 - N1v) * N2 + k2 + i].imag();
			}
		}
		if (N1 <= SFFT_NMAX) {
			sfft_complex((fft_simd4*) zv.data(), N1);
		} else {
			fft_scramble(zv.data(), N1);
			fft(zv.data(), N1);
		}
		for (int k1 = 0; k1 < N1; k1++) {
			for (int i = 0; i < SIMD_SIZE; i++) {
				X[k1 * N2 + k2 + i].real() = zv[k1].real()[i];
				X[k1 * N2 + k2 + i].imag() = zv[k1].imag()[i];
			}
		}
	}
	for (int k2 = N2v; k2 < N2; k2++) {
		for (int n1r = 0; n1r < N1voS; n1r++) {
			for (int n1c = 0; n1c < SIMD_SIZE; n1c++) {
				const int n1 = n1r * SIMD_SIZE + n1c;
				zs[n1].real() = Yv[n1r * N2 + k2].real()[n1c];
				zs[n1].imag() = Yv[n1r * N2 + k2].imag()[n1c];
			}
		}
		for (int n1 = N1v; n1 < N1; n1++) {
			zs[n1] = Ys[(n1 - N1v) * N2 + k2];
		}
		if (N1 <= SFFT_NMAX) {
			sfft_complex((double*) zs.data(), N1);
		} else {
			fft(zs.data(), N1);
		}
		for (int k1 = 0; k1 < N1; k1++) {
			X[k1 * N2 + k2] = zs[k1];
		}
	}
	sws.destroy(std::move(zs));
	sws.destroy(std::move(Ys));
	vws.destroy(std::move(Yv));
	vws.destroy(std::move(zv));
}

std::vector<fft_method> possible_start_ffts(int N) {
	std::vector<fft_method> ffts;
	fft_method m;
	if (is_prime(N) && N > SFFT_NMAX) {
		m.type = FFT_RADERS_PADDED;
		ffts.push_back(m);
		m.type = FFT_RADERS;
		ffts.push_back(m);
	}
	bool found = false;
	for (int R = SIMD_SIZE; R <= SFFT_NMAX; R++) {
		if (N % R == 0) {
			m.type = FFT_CT;
			m.R = R;
			ffts.push_back(m);
			found = true;
		}
	}
	if (!found) {
		m.R = SFFT_NMAX + 1;
		while (N % m.R != 0 && m.R < N) {
			m.R++;
		}
		if (m.R < N) {
			m.type = FFT_CT;
			ffts.push_back(m);
		}
	}
	if (!padded_length(N) && (is_prime(N) || (N % 2 == 0 && is_prime(N / 2)) || (N % 3 == 0 && is_prime(N / 3)))) {
		m.type = FFT_BLUE;
		ffts.push_back(m);
	}
	return ffts;
}

std::vector<fft_method> possible_ffts(int N) {
	std::vector<fft_method> ffts;
	fft_method m;
	if (lround(sqrt(N)) * lround(sqrt(N)) == N) {
		m.type = FFT_6;
		ffts.push_back(m);
		return ffts;
	}
	if (N % 2 == 0) {
		for (m.R = 4; m.R <= std::min(N, SFFT_NMAX); m.R += 4) {
			if (N % m.R == 0) {
				m.type = FFT_SPLIT;
				ffts.push_back(m);
			}
		}
		for (m.R = 4; m.R <= std::min(N, SFFT_NMAX); m.R += 4) {
			if (N % m.R == 0) {
				m.type = FFT_SPLIT_CONJ;
				ffts.push_back(m);
			}
		}
	}
	for (m.R = 2; m.R <= std::min(N, SFFT_NMAX); m.R++) {
		if (N % m.R == 0) {
			m.type = FFT_CT;
			ffts.push_back(m);
			m.type = FFT_CONJ;
			ffts.push_back(m);
		}
	}
	if (ffts.size() == 0) {
		auto factors = prime_factorization(N);
		if (factors.size() == 1) {
			m.type = FFT_RADERS;
			m.R = N;
			ffts.push_back(m);
			m.type = FFT_RADERS_PADDED;
			ffts.push_back(m);
		} else {
			for (auto i = factors.begin(); i != factors.end(); i++) {
				m.R = i->first;
				m.type = FFT_CT;
				ffts.push_back(m);
				m.type = FFT_CONJ;
				ffts.push_back(m);
			}
		}
	}
	return ffts;
}

void fft(const fft_method& method, complex<fft_simd4>* X, int N) {
	switch (method.type) {
	case FFT_SPLIT:
		return fft_split(method.R, X, N);
	case FFT_SPLIT_CONJ:
		return fft_split_conjugate(method.R, X, N);
	case FFT_CT:
		return fft_cooley_tukey(method.R, X, N);
	case FFT_CONJ:
		return fft_conjugate(method.R, X, N);
	case FFT_6:
		return fft_six_step(X, N);
	case FFT_RADERS:
		return fft_raders(X, N);
	case FFT_RADERS_PADDED:
		return fft_raders_padded(X, N);
	}
}

void fft(const fft_method& method, complex<double>* X, int N) {
	switch (method.type) {
	case FFT_CT:
		return fft2simd(method.R, X, N);
	case FFT_BLUE:
		return fft_bluestein(X, N);
	case FFT_RADERS:
		return fft_raders(X, N, false);
	case FFT_RADERS_PADDED:
		return fft_raders(X, N, true);
	}
}

void fft_indices(const fft_method& method, int* I, int N) {
	switch (method.type) {
	case FFT_SPLIT:
		return fft_split_indices(method.R, I, N);
	case FFT_SPLIT_CONJ:
		return fft_split_conjugate_indices(method.R, I, N);
	case FFT_CT:
		return fft_cooley_tukey_indices(method.R, I, N);
	case FFT_CONJ:
		return fft_conjugate_indices(method.R, I, N);
	case FFT_RADERS:
		return fft_raders_indices(I, N);
	case FFT_RADERS_PADDED:
		return fft_raders_padded_indices(I, N);
	case FFT_6:
		return fft_six_step_indices(I, N);
	}
}

fft_method select_fft(int N) {
	static std::unordered_map<int, fft_method> cache;
	auto iter = cache.find(N);
	if (iter == cache.end()) {
		std::vector<complex<fft_simd4>> X(N);
		const auto tests = possible_ffts(N);
		const int M = tests.size();
		std::vector<double> timers(M);
		for (int m = 0; m < M; m++) {
			for (int n = 0; n < N; n++) {
				X[n].real() = rand1();
				X[n].imag() = rand1();
			}
			fft(tests[m], X.data(), N);
			std::vector<double> times;
			for (int n = 0; n < NTRIAL; n++) {
				timer tm;
				for (int n = 0; n < N; n++) {
					X[n].real() = rand1();
					X[n].imag() = rand1();
				}
				tm.start();
				fft(tests[m], X.data(), N);
				tm.stop();
				times.push_back(tm.read());
			}
			std::sort(times.begin(), times.end());
			timers[m] = times[(times.size() + 1) / 2];
			//.		printf("%i - %s - %e\n", N, fft_method_string(tests[m]).c_str(), timers[m]);
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

void fft(complex<fft_simd4>* X, int N) {
	if (N <= SFFT_NMAX) {
		sfft_complex((fft_simd4*) X, N);
		return;
	}
	fft(select_fft(N), X, N);
}

fft_method select_start_fft(int N) {
	static std::unordered_map<int, fft_method> cache;
	auto iter = cache.find(N);
	if (iter == cache.end()) {
		std::vector<complex<double>> X(N);
		const auto tests = possible_start_ffts(N);
		const int M = tests.size();
		std::vector<double> timers(M);
		for (int m = 0; m < M; m++) {
			for (int n = 0; n < N; n++) {
				X[n].real() = rand1();
				X[n].imag() = rand1();
			}
			fft(tests[m], X.data(), N);
			std::vector<double> times;
			for (int n = 0; n < NTRIAL; n++) {
				timer tm;
				for (int n = 0; n < N; n++) {
					X[n].real() = rand1();
					X[n].imag() = rand1();
				}
				tm.start();
				fft(tests[m], X.data(), N);
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

void fft(complex<double>* X, int N) {
	if (N <= SFFT_NMAX) {
		sfft_complex((double*) X, N);
		return;
	}
	fft(select_start_fft(N), X, N);
}

void fft_indices(int* I, int N) {
	if (N <= SFFT_NMAX) {
		return;
	}
	fft_indices(select_fft(N), I, N);
}

const std::vector<int>& fft_indices(int N) {
	static std::unordered_map<int, std::shared_ptr<std::vector<int>>>cache;
	auto iter = cache.find(N);
	if (iter == cache.end()) {
		std::vector<int> I(N);
		std::vector<int> J(N);
		std::iota(I.begin(), I.end(), 0);
		fft_indices(I.data(), N);
		cache[N] = std::make_shared<std::vector<int>>(std::move(I));
		iter = cache.find(N);
	}
	return *iter->second;
}

const std::vector<int>& fft_inv_indices(int N) {
	static std::unordered_map<int, std::shared_ptr<std::vector<int>>>cache;
	auto iter = cache.find(N);
	if (iter == cache.end()) {
		std::vector<int> I(N);
		std::vector<int> J(N);
		std::iota(I.begin(), I.end(), 0);
		fft_indices(I.data(), N);
		for (int n = 0; n < N; n++) {
			J[I[n]] = n;
		}
		cache[N] = std::make_shared<std::vector<int>>(std::move(J));
		iter = cache.find(N);
	}
	return *iter->second;
}

template<class T>
void fft_scramble1(complex<T>* X, int N) {
	if (N <= SFFT_NMAX) {
		return;
	}
	static std::unordered_map<int, std::shared_ptr<std::vector<int>>>cache;
	auto iter = cache.find(N);
	if (iter == cache.end()) {
		std::vector<int> I = fft_inv_indices(N);
		cache[N] = std::make_shared<std::vector<int>>(std::move(I));
		iter = cache.find(N);
	}
	const auto& I = *iter->second;
	fft_permute(I, X);
}

template<class T>
void fft_scramble2(complex<T>* X, int N) {
	if (N <= SFFT_NMAX) {
		return;
	}
	static std::unordered_map<int, std::shared_ptr<std::vector<int>>>cache;
	auto iter = cache.find(N);
	if (iter == cache.end()) {
		std::vector<int> I = fft_inv_indices(N);
		for (int n = 1; n < N - n; n++) {
			std::swap(I[n], I[N - n]);
		}
		cache[N] = std::make_shared<std::vector<int>>(std::move(I));
		iter = cache.find(N);
	}
	const auto& I = *iter->second;
	fft_permute(I, X);
}

void fft_scramble(complex<double>* X, int N) {
	fft_scramble1(X, N);
}

void fft_scramble(complex<fft_simd4>* X, int N) {
	fft_scramble1(X, N);
}

void fft_scramble_inv(complex<double>* X, int N) {
	fft_scramble2(X, N);
}

void fft_scramble_inv(complex<fft_simd4>* X, int N) {
	fft_scramble2(X, N);
}

