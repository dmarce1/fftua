#include <vector>
#include <numeric>
#include <cstdio>
#include <cstdlib>
#include <unordered_map>
#include <map>
#include <cmath>
#include <set>
#include <complex>
#include <fftw3.h>
#include <cassert>
#include <memory>
#include <stack>
#include "types.hpp"
#include "util.hpp"

std::map<int, int> prime_factorization(int N);

bool are_coprime(int a, int b);

class indexer {
	std::vector<int> N;
	std::vector<int> Nprod;
public:
	indexer(std::vector<int>&& n) {
		N = std::move(n);
		Nprod.resize(N.size());
		Nprod[0] = 1;
		for (unsigned n = 0; n < N.size() - 1; n++) {
			Nprod[n + 1] = Nprod[n] * N[n];
		}
	}
	int index(int index, int dim) const {
		return (index / Nprod[dim]) % N[dim];
	}
	int max() const {
		return Nprod.back() * N.back();
	}
};

int mod(int a, int b) {
	while (a < 0) {
		a += b;
	}
	return a % b;
}

bool power_of(int N, int M) {
	while (N > 1) {
		if (N % M != 0) {
			return false;
		}
		N /= M;
	}
	return true;
}

int mod_pow(int a, int b, int m) {
	int rc = 1;
	int apow = a;
	if (a < 0 && b % 2 == 1) {
		return -mod_pow(-a, b, m);
	}
	while (b) {
		if (b & 1) {
			rc = ((long long) (rc % m) * (long long) (apow % m)) % m;
		}
		b >>= 1;
		apow = ((long long) (apow % m) * (long long) (apow % m)) % m;
	}
	return rc;
}

int generator(int N) {
	static thread_local std::unordered_map<int, int> values;
	auto i = values.find(N);
	auto factors = prime_factorization(N);
	int P = factors.begin()->first;
	int c = factors.begin()->second;
	int M = std::pow(P, c - 1) * (P - 1);
	if (i == values.end()) {
		for (int g = 2;; g++) {
			std::set<int> I;
			bool fail = false;
			for (int m = 1; m < M; m++) {
				if (are_coprime(m, N)) {
					int n = mod_pow(g, m, N);
					if (I.find(n) == I.end()) {
						I.insert(n);
					} else {
						fail = true;
						break;
					}
				}
			}
			if (!fail) {
				values[N] = g;
				i = values.find(N);
				break;
			}
		}
	}
	return i->second;
}

std::vector<int> raders_ginvq(int N) {
	static thread_local std::unordered_map<int, std::shared_ptr<std::vector<int>>>values;
	auto i = values.find(N);
	if (i == values.end()) {
		const int g = generator(N);
		std::vector<int> ginvq;
		for (int q = 0; q < N - 1; q++) {
			ginvq.push_back(mod_inv(mod_pow(g, q, N), N));
		}
		ginvq.push_back(0);
		values[N] = std::make_shared<std::vector<int>>(std::move(ginvq));
		i = values.find(N);
	}
	return *i->second;
}

const std::vector<int> raders_gq(int N) {
	static thread_local std::unordered_map<int, std::shared_ptr<std::vector<int>>>values;
	auto i = values.find(N);
	if (i == values.end()) {
		const int g = generator(N);
		std::vector<int> gq;
		for (int q = 0; q < N - 1; q++) {
			gq.push_back(mod_pow(g, q, N));
		}
		gq.push_back(0);
		values[N] = std::make_shared<std::vector<int>>(std::move(gq));
		i = values.find(N);
	}
	return *i->second;
}

double fftw(std::vector<complex<double>>& x) {
	const int N = x.size();
	static std::unordered_map<int, fftw_plan> plans;
	static std::unordered_map<int, fftw_complex*> in;
	static std::unordered_map<int, fftw_complex*> out;
	if (plans.find(N) == plans.end()) {
		in[N] = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
		out[N] = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
		//out[N] = in[N];
		plans[N] = fftw_plan_dft_1d(N, in[N], out[N], FFTW_FORWARD, FFTW_MEASURE | FFTW_NO_SIMD);
	}
	auto* i = in[N];
	auto* o = out[N];
	for (int n = 0; n < N; n++) {
		i[n][0] = x[n].real();
		i[n][1] = x[n].imag();
	}
	timer tm;
	tm.start();
	fftw_execute(plans[N]);
	tm.stop();
	for (int n = 0; n < N; n++) {
		x[n].real() = (o[n][0]);
		x[n].imag() = (o[n][1]);
	}
	return tm.read();
}

const std::vector<complex<double>>& twiddles(int N) {
	using entry_type = std::shared_ptr<std::vector<complex<double>>>;
	static std::unordered_map<int, entry_type> cache;
	auto iter = cache.find(N);
	if (iter != cache.end()) {
		return *(iter->second);
	} else {
		std::vector<complex<double>> W(round_up(N, SIMD_SIZE));
		for (int n = 0; n < N; n++) {
			W[n].real() = cos(-2.0 * M_PI * n / N);
			W[n].imag() = sin(-2.0 * M_PI * n / N);
		}
		cache[N] = std::make_shared<std::vector<complex<double>>>(std::move(W));
		return *(cache[N]);
	}
}

const std::vector<std::vector<double>>& cos_twiddles(int N1, int N2) {
	using entry_type = std::shared_ptr<std::vector<std::vector<double>>>;
	static std::unordered_map<int, std::unordered_map<int, entry_type>> cache;
	const int N = N1 * N2;
	auto iter = cache[N1].find(N2);
	if (iter != cache[N1].end()) {
		return *(iter->second);
	} else {
		std::vector<std::vector<double>> W(N1, std::vector<double>(N2));
		for (int n1 = 0; n1 < N1; n1++) {
			for (int n2 = 0; n2 < N2; n2++) {
				W[n1][n2] = cos(2.0 * M_PI * n1 * n2 / N);
			}
		}
		cache[N1][N2] = std::make_shared<std::vector<std::vector<double>>>(std::move(W));
		return *(cache[N1][N2]);
	}
}

const std::vector<std::vector<double>>& sin_twiddles(int N1, int N2) {
	using entry_type = std::shared_ptr<std::vector<std::vector<double>>>;
	static std::unordered_map<int, std::unordered_map<int, entry_type>> cache;
	const int N = N1 * N2;
	auto iter = cache[N1].find(N2);
	if (iter != cache[N1].end()) {
		return *(iter->second);
	} else {
		std::vector<std::vector<double>> W(N1, std::vector<double>(N2));
		for (int n1 = 0; n1 < N1; n1++) {
			for (int n2 = 0; n2 < N2; n2++) {
				W[n1][n2] = sin(-2.0 * M_PI * n1 * n2 / N);
			}
		}
		cache[N1][N2] = std::make_shared<std::vector<std::vector<double>>>(std::move(W));
		return *(cache[N1][N2]);
	}
}

const std::vector<std::vector<complex<fft_simd4>>>& vector_twiddles(int N1, int N2) {
	using entry_type = std::shared_ptr<std::vector<std::vector<complex<fft_simd4>>>>;
	static std::unordered_map<int, std::unordered_map<int, entry_type>> cache;
	auto iter = cache[N1].find(N2);
	if (iter != cache[N1].end()) {
		return *(iter->second);
	} else {
		const int N = N1 * N2;
		const int N1v = round_down(N1, SIMD_SIZE);
		std::vector<std::vector<complex<fft_simd4>>>W(N1v/SIMD_SIZE, std::vector<complex<fft_simd4>>(N2));
		for (int n = 0; n < N1v; n += SIMD_SIZE) {
			for (int k = 0; k < N2; k++) {
				for (int i = 0; i < SIMD_SIZE; i++) {
					W[n / SIMD_SIZE][k].real()[i] = cos(-2.0 * M_PI * (n+i) * k / N);
					W[n / SIMD_SIZE][k].imag()[i] = sin(-2.0 * M_PI * (n+i) * k / N);
				}
			}
		}
		cache[N1][N2] = std::make_shared<std::vector<std::vector<complex<fft_simd4>>>>(std::move(W));
		return *(cache[N1][N2]);
	}
}

const std::vector<complex<fft_simd4>>& vector_twiddles2(int N1, int N2) {
	using entry_type = std::shared_ptr<std::vector<complex<fft_simd4>>>;
	static std::unordered_map<int, std::unordered_map<int, entry_type>> cache;
	auto iter = cache[N1].find(N2);
	if (iter != cache[N1].end()) {
		return *(iter->second);
	} else {
		const int N = N1 * N2;
		std::vector<complex<fft_simd4>> W(N1 * N2 / SIMD_SIZE);
		for (int n = 0; n < N1; n += SIMD_SIZE) {
			for (int k = 0; k < N2; k++) {
				for (int i = 0; i < SIMD_SIZE; i++) {
					W[(k * N1 + n) / SIMD_SIZE].real()[i] = cos(-2.0 * M_PI * (n+i) * k / N);
					W[(k * N1 + n) / SIMD_SIZE].imag()[i] = sin(-2.0 * M_PI * (n+i) * k / N);
				}
			}
		}
		cache[N1][N2] = std::make_shared<std::vector<complex<fft_simd4>>>(std::move(W));
		return *(cache[N1][N2]);
	}
}

const std::vector<std::vector<complex<fft_simd4>>>& shifted_vector_twiddles(int N1, int N2) {
	std::swap(N1,N2);
	using entry_type = std::shared_ptr<std::vector<std::vector<complex<fft_simd4>>>>;
	static std::unordered_map<int, std::unordered_map<int, entry_type>> cache;
	auto iter = cache[N1].find(N2);
	if (iter != cache[N1].end()) {
		return *(iter->second);
	} else {
		const int N = N1 * N2;
		const int N1v = round_down(N1, SIMD_SIZE);
		std::vector<std::vector<complex<fft_simd4>>>W(N1v/SIMD_SIZE, std::vector<complex<fft_simd4>>(N2));
		for (int n = 0; n < N1v; n += SIMD_SIZE) {
			for (int k = 0; k < N2; k++) {
				for (int i = 0; i < SIMD_SIZE; i++) {
					W[n / SIMD_SIZE][k].real()[i] = cos(-2.0 * M_PI * (n+i+1) * k / N);
					W[n / SIMD_SIZE][k].imag()[i] = sin(-2.0 * M_PI * (n+i+1) * k / N);
				}
			}
		}
		cache[N1][N2] = std::make_shared<std::vector<std::vector<complex<fft_simd4>>>>(std::move(W));
		return *(cache[N1][N2]);
	}
}

bool padded_length(int M) {
	auto factors = prime_factorization(M);
	const int two_pow = factors.begin()->first == 2 ? factors.begin()->second : 0;
	for (auto i = factors.begin(); i != factors.end(); i++) {
		if ((i->second > two_pow) || (i->first != 2 && i->first != 3 && i->first != 5)) {
			return false;
		}
	}
	return true;
}

int compute_padding(int N) {
	static std::unordered_map<int, int> cache;
	auto iter = cache.find(N);
	if (iter == cache.end()) {
		int M = 2 * N - 1;
		bool done;
		do {
			if (!padded_length(M)) {
				M++;
				done = false;
			} else {
				done = true;
			}
		} while (!done);
		cache[N] = M;
		iter = cache.find(N);
	}
	return iter->second;
}

const std::vector<complex<double>>& raders_twiddle(int N, int M) {
	using entry_type = std::shared_ptr<std::vector<complex<double>>>;
	static std::unordered_map<int, std::unordered_map<int, entry_type>> cache;
	auto iter = cache[N].find(M);
	if (iter != cache[N].end()) {
		return *(iter->second);
	} else {
		auto factors = prime_factorization(N);
		assert(factors.size() == 1);
		int P = factors.begin()->first;
		int c = factors.begin()->second;
		int L = std::pow(P, c - 1) * (P - 1);
		std::vector<complex<double>> b(M, 0.0);
		const auto tws = twiddles(N);
		const auto ginvq = raders_ginvq(N);
		for (int q = 0; q < L; q++) {
			b[q] = (1.0 / M) * tws[ginvq[q]];
		}
		if (M != N - 1) {
			for (int q = 1; q < L; q++) {
				b[M - q] = b[L - q];
			}
		}
		fftw(b);
		b.resize(round_up(b.size(), SIMD_SIZE));
		cache[N][M] = std::make_shared<std::vector<complex<double>>>(std::move(b));
		return *(cache[N][M]);
	}
}

void fftw_dht(std::vector<double>& xin);

const std::vector<double>& raders_twiddle_real(int N, int M) {
	using entry_type = std::shared_ptr<std::vector<double>>;
	static std::unordered_map<int, std::unordered_map<int, entry_type>> cache;
	auto iter = cache[N].find(M);
	if (iter != cache[N].end()) {
		return *(iter->second);
	} else {
		auto factors = prime_factorization(N);
		assert(factors.size() == 1);
		int P = factors.begin()->first;
		int c = factors.begin()->second;
		int L = std::pow(P, c - 1) * (P - 1);
		std::vector<double> b(M, 0.0);
		const auto tws = twiddles(N);
		const auto ginvq = raders_ginvq(N);
		for (int q = 0; q < L; q++) {
			b[q] = (1.0 / M) * (tws[ginvq[q]].real() - tws[ginvq[q]].imag());
		}
		if (M != N - 1) {
			for (int q = 1; q < L; q++) {
				b[M - q] = b[L - q];
			}
		}
		fftw_dht(b);
		b.resize(round_up(b.size(), SIMD_SIZE));
		cache[N][M] = std::make_shared<std::vector<double>>(std::move(b));
		return *(cache[N][M]);
	}
}

const std::vector<complex<double>>& bluestein_multiplier(int N) {
	static thread_local std::unordered_map<int, std::shared_ptr<std::vector<complex<double>>> >values;
	auto i = values.find(N);
	if (i == values.end()) {
		std::vector<complex<double>> z(N);
		for (int n = 0; n < N; n++) {
			const auto phi = fmodl(M_PIl * n * n / N, 2.0l * M_PIl);
			z[n].real() = cos(phi);
			z[n].imag() = -sin(phi);
		}
		values[N] = std::make_shared<std::vector<complex<double>>>(std::move(z));
		i = values.find(N);
	}
	return *i->second;
}

const std::vector<complex<double>>& bluestein_filter(int N, int M) {
	static thread_local std::unordered_map<int, std::unordered_map<int, std::shared_ptr<std::vector<complex<double>>>>>values;
	auto i = values[N].find(M);
	if (i == values[N].end()) {
		std::vector<complex<double>> z(M, 0.0);
		for (int n = 0; n < N; n++) {
			const auto phi = fmodl(M_PIl * n * n / N, 2.0l * M_PIl);
			z[n].real() = cos(phi) / M;
			z[n].imag() = sin(phi) / M;
		}
		for (int n = 1; n < N; n++) {
			z[M - n] = z[n];
		}
		fftw(z);
		values[N][M] = std::make_shared<std::vector<complex<double>>>(std::move(z));
		i = values[N].find(M);
	}
	return *i->second;
}

__int128 factorial(__int128 k) {
	if (k <= 1) {
		return 1;
	} else {
		return k * factorial(k - 1);
	}
}

std::vector<std::vector<int>> nchoosek(int n, int k) {
	std::vector<std::vector<int>> rc;
	std::vector<int> combo(k);
	std::iota(combo.begin(), combo.end(), 0);
	if (n == k) {
		rc.push_back(combo);
	} else {
		bool done = false;
		while (!done) {
			rc.push_back(combo);
			int dim = k - 1;
			while (combo[dim] == dim + n - k) {
				dim--;
				if (dim < 0) {
					done = true;
					break;
				}
			}
			if (!done) {
				combo[dim]++;
				for (int i = dim + 1; i < k; i++) {
					combo[i] = combo[i - 1] + 1;
				}
			}
		}
	}
	rc.push_back(combo);
	return rc;
}

bool is_prime(int n) {
	static thread_local std::unordered_map<int, bool> values;
	auto i = values.find(n);
	if (i == values.end()) {
		bool v = true;
		if (n == 1) {
			v = false;
		} else {
			int kmax = sqrt(n);
			for (auto i = 2; i <= kmax; i++) {
				if (n % i == 0) {
					v = false;
					break;
				}
				if (i > kmax) {
					break;
				}
			}
		}
		values[n] = v;
		i = values.find(n);
	}
	return i->second;
}

int greatest_prime_factor(int N) {
	static thread_local std::unordered_map<int, int> values;
	auto i = values.find(N);
	if (i == values.end()) {
		int v;
		for (int n = 2; n <= N; n++) {
			if (N % n == 0 && is_prime(n)) {
				v = n;
			}
		}
		values[N] = v;
		i = values.find(N);
	}
	return i->second;
}

int least_prime_factor(int N) {
	static thread_local std::unordered_map<int, int> values;
	auto i = values.find(N);
	if (i == values.end()) {
		int v;
		for (int n = N; n >= 2; n--) {
			if (N % n == 0 && is_prime(n)) {
				v = n;
			}
		}
		values[N] = v;
		i = values.find(N);
	}
	return i->second;
}

std::map<int, int> prime_factorization(int N) {
	static thread_local std::unordered_map<int, std::map<int, int>> values;
	auto i = values.find(N);
	if (i == values.end()) {
		std::map<int, int> map;
		while (N != 1) {
			int k = greatest_prime_factor(N);
			if (map.find(k) == map.end()) {
				map[k] = 0;
			}
			map[k]++;
			N /= k;
		}
		values[N] = std::move(map);
		i = values.find(N);
	}
	return i->second;
}

int totient(int N) {
	auto P = prime_factorization(N);
	int T = 1;
	for (auto i = P.begin(); i != P.end(); i++) {
		T *= std::pow(i->first, i->second - 1);
		T *= (i->first - 1);
	}
	return T;
}

int mod_inv(int a, int m) {
	return mod_pow(a, totient(m) - 1, m);
}

bool are_coprime(int a, int b) {
	auto afacs = prime_factorization(a);
	auto bfacs = prime_factorization(b);
	for (auto i : afacs) {
		for (auto j : bfacs) {
			if (i.first == j.first) {
				return false;
			}
		}
	}
	return true;
}

double fftw_real(std::vector<complex<double>>& xout, const std::vector<double>& xin) {
	const int N = xin.size();
	static std::unordered_map<int, fftw_plan> plans;
	static std::unordered_map<int, double*> in;
	static std::unordered_map<int, fftw_complex*> out;
	if (plans.find(N) == plans.end()) {
		in[N] = (double*) malloc(sizeof(double) * N);
		out[N] = (fftw_complex*) malloc(sizeof(fftw_complex) * (N / 2 + 1));
		plans[N] = fftw_plan_dft_r2c_1d(N, in[N], out[N], FFTW_MEASURE );
	}
	auto* i = in[N];
	auto* o = out[N];
	for (int n = 0; n < N; n++) {
		i[n] = xin[n];
	}
	timer tm;
	tm.start();
	fftw_execute(plans[N]);
	tm.stop();
	for (int n = 0; n < N / 2 + 1; n++) {
		xout[n].real() = (o[n][0]);
		xout[n].imag() = (o[n][1]);
	}
	return tm.read();
}

void fftw_dht(std::vector<double>& h) {
	int N = h.size();
	std::vector<complex<double>> Y(N / 2 + 1);
	fftw_real(Y, h);
	for (int n = 1; n < N - n; n++) {
		h[n] = Y[n].real() - Y[n].imag();
		h[N - n] = Y[n].real() + Y[n].imag();
	}
	h[0] = Y[0].real();
	if (N % 2 == 0) {
		h[N / 2] = Y[N / 2].real();
	}
}

