#include "util.hpp"
#include "timer.hpp"
#include <fftw3.h>
#include <cmath>
#include <memory>
#include <unordered_map>

double fftw(std::vector<complex<double>>& x) {
	timer tm;
	const int N = x.size();
	static std::unordered_map<int, fftw_plan> plans;
	static std::unordered_map<int, fftw_complex*> in;
	static std::unordered_map<int, fftw_complex*> out;
	if (plans.find(N) == plans.end()) {
		in[N] = (fftw_complex*) malloc(sizeof(fftw_complex) * N);
		out[N] = (fftw_complex*) malloc(sizeof(fftw_complex) * N);
		plans[N] = fftw_plan_dft_1d(N, in[N], out[N], FFTW_FORWARD, FFTW_MEASURE | FFTW_NO_SIMD);
	}
	auto* i = in[N];
	auto* o = out[N];
	for (int n = 0; n < N; n++) {
		i[n][0] = x[n].real();
		i[n][1] = x[n].imag();
	}
	tm.start();
	fftw_execute(plans[N]);
	tm.stop();
	for (int n = 0; n < N; n++) {
		x[n].real() = o[n][0];
		x[n].imag() = o[n][1];
	}
	return tm.read();
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
		int v = -1;
		int n = N;
		while (true) {
			if (N % n == 0 && is_prime(n)) {
				v = n;
				break;
			}
			n--;
		}
		values[N] = v;
		i = values.find(N);
	}
	return i->second;
}

const std::map<int, int>& prime_factorization(int N) {
	static thread_local std::unordered_map<int, std::shared_ptr<std::map<int, int>>>values;
	auto i = values.find(N);
	if (i == values.end()) {
		std::map<int, int> map;
		int n = N;
		if( n % 2 == 0 ) {
			map[2] = 0;
			while( n % 2 == 0 ) {
				map[2]++;
				n >>= 1;
			}
		}
		int r = 3;
		while( r <= sqrt(N)) {
			if( n % r == 0 ) {
				map[r] = 0;
				while( n % r == 0 ) {
					map[r]++;
					n /= r;
				}
			}
			r += 2;
		}
		if( n > 2 ) {
			map[N] = 1;
		}
		values[N] = std::make_shared<std::map<int, int>>(std::move(map));
		i = values.find(N);
	}
	return *i->second;
}

const std::vector<complex<double>>& twiddles(int N) {
	static thread_local std::unordered_map<int, std::shared_ptr<std::vector<complex<double>>> >values;
	auto i = values.find(N);
	if (i == values.end()) {
		std::vector<complex<double>> W(N);
		for( int n = 0; n < N; n++) {
			const double theta = -2.0 * M_PI * n / N;
			W[n].real() = cos(theta);
			W[n].imag() = sin(theta);
		}
		values[N] = std::make_shared<std::vector<complex<double>>>(std::move(W));
		i = values.find(N);
	}
	return *i->second;
}

double rand1() {
	return (rand() + 0.5) / (RAND_MAX + 1.0);
}
