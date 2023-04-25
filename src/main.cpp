#include <complex>
#include <vector>
#include <unordered_map>
#include <memory>
#include <chrono>
#include <numeric>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include "sfft.hpp"
#include "util.hpp"
#include "types.hpp"
#include "fftu.hpp"

void * operator new(std::size_t n) {
	void* memptr;
	posix_memalign(&memptr, 32, round_up(n, 32));
	return memptr;
}

void operator delete(void * p) {
	free(p);
}

void *operator new[](std::size_t n) {
	void* memptr;
	posix_memalign(&memptr, 32, round_up(n, 32));
	return memptr;
}

void operator delete[](void *p) {
	free(p);
}

#define SIMD_SIZE 4

double FFT(std::vector<complex<double>>& X) {
	int N = X.size();
	static std::vector<complex<double>> Y;
	Y.resize(N);
	timer tm;
	tm.start();
	fft(X.data(), N);
	tm.stop();
	return tm.read();
}

int permute_index(int index, int width) {
	for (int pos = 0; pos < width; pos++) {
		if (pos < width - 1) {
			int lo = (1 << pos) - 1;
			int hi = (~lo) << 2;
			int bits = (index & (0x3 << pos)) >> pos;
			index = (index & hi) | ((index & lo) << 2) | bits;
			pos++;
		} else {
			int lo = (1 << pos) - 1;
			int hi = (~lo) << 1;
			int bit = (index & (1 << pos)) >> pos;
			index = (index & hi) | ((index & lo) << 1) | bit;
		}
	}
	return index;
}

#include <fenv.h>
void fft_scramble_real(double* X, int N);

int main(int argc, char **argv) {
	//printf( "PRIMITIVE ROOT OF 93871 = %i\n", generator(93871));
//	feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
	timer tm3, tm4;
	constexpr int N = 16;
	std::vector<complex<double>> x(N);
	std::vector<complex<double>> Y(N / 2 + 1);
	std::vector<double> y(N);
	for (int n = 0; n < N; n++) {
		x[n].real() = rand1();
		x[n].imag() = rand1();
	}
	x[0].imag() = x[N / 2].imag() = 0.0;
	for (int n = 1; n < N - n; n++) {
		if (n % 2 == 0) {
			x[N - n] = x[n].conj();
		} else {
			x[N - n] = -x[n].conj();
		}
	}
	y[0] = x[0].real();
	if (N % 2 == 0) {
		y[N / 2] = x[N / 2].real();
	}
	for (int n = 1; n < N / 2; n++) {
		auto e = y[2 * n];
		auto o = y[2 * n + 1];
		y[n] = e;
		if( n != 0 ) {
			y[N - n] = e;
		}
	}
	fftw_real(Y, y);
	fftw(x);
	auto Z = Y;
	for (int n = 0; n < N - n; n++) {
		Y[n].real() = (Z[n].real() + Z[n].imag());
		Y[n].imag() = (Z[n].real() - Z[n].imag());
	}
	for (int n = 0; n < N - n; n++) {
		printf("%i %e %e | %e %e\n", n, Y[n].real(), Y[n].imag(), x[n].real(), x[n].imag());
	}
		double t3 = 0.0;
	double t4 = 0.0;
	std::vector<int> Ns;
	double score = 0.0;
	int cnt = 0;
	for (int N = 13; N <= 1024 * 1024; N = 11 * N / 10) {
		auto pfac = prime_factorization(N);
		if (pfac.size() != 1) {
			continue;
		}
		double avg_err = 0.0;
		double t1 = 0.0;
		double t2 = 0.0;

		for (int i = 0; i < 2; i++) {
			std::vector<double> x(N);
			std::vector<double> y(N);
			std::vector<complex<double>> X(N / 2 + 1);
			std::vector<complex<double>> Y(N / 2 + 1);
			for (int n = 0; n < N; n++) {
				x[n] = (y[n] = rand1());
			}
			x[1] = y[1] = 1.0;
			int N1 = std::pow(pfac.rbegin()->first, pfac.rbegin()->second);
			if (i == 0) {
				fftw_real(Y, y);
				fft_raders_real(x.data(), N, false);
			} else {
				auto b = fftw_real(Y, y);
				timer tm;
				tm.start();
				//		fft_scramble_real(x.data(), N);
				fft_raders_real(x.data(), N, false);
				X[0].real() = x[0];
				X[0].imag() = 0.0;
				for (int n = 1; n < N - n; n++) {
					X[n].real() = x[n];
					X[n].imag() = x[N - n];
				}
				if (N % 2 == 0) {
					X[N / 2].real() = x[N / 2];
					X[N / 2].imag() = 0.0;
				}
				tm.stop();
				t2 += tm.read();
				t4 += tm.read();
				t1 += b;
				t3 += b;
				for (int n = 0; n < N / 2 + 1; n++) {
					double x = X[n].real() - Y[n].real();
					double y = X[n].imag() - Y[n].imag();
					double err = sqrt(x * x + y * y);
					avg_err += err;
					printf("%16e %16e | %16e %16e | %16e %16e\n", X[n].real(), X[n].imag(), Y[n].real(), Y[n].imag(), X[n].real() - Y[n].real(), X[n].imag() - Y[n].imag());
				}
			}
		}
		abort();
		std::string f;
		for (auto i = pfac.begin(); i != pfac.end(); i++) {
			f += "(" + std::to_string(i->first) + "^" + std::to_string(i->second) + ")";
		}
		printf("%i: %32s ", N, f.c_str());
		fflush(stdout);
		avg_err /= (20 * N);
		score *= cnt;
		score += t1 / (t2 + 1e-20);
		cnt++;
		score /= cnt;
		printf("%c| %e %e %e %e %e | %e\n", (pfac.size() == 1 && pfac.begin()->second == 1) ? '*' : ' ', avg_err, t1, t2, t1 / (t2 + 1e-20), t4, score);
	}
	return 0;
}
