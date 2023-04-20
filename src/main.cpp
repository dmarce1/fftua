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

int main(int argc, char **argv) {
	//printf( "PRIMITIVE ROOT OF 93871 = %i\n", generator(93871));
//	feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
	timer tm3, tm4;
	double t3 = 0.0;
	double t4 = 0.0;
	std::vector<int> Ns;
	double score = 0.0;
	int cnt = 0;
	for (int N = 5; N <= 16 * 1024 * 1024; N *= 5) {
		double avg_err = 0.0;
		double t1 = 0.0;
		double t2 = 0.0;
		//	for (int i = 0; i < 21; i++) {
		std::vector<double> x(N);
		std::vector<double> y(N);
		std::vector<complex<double>> X(N / 2 + 1);
		std::vector<complex<double>> Y(N / 2 + 1);
		for (int n = 0; n < N; n++) {
			x[n] = (y[n] = rand1());
		}
		/*		if (i == 0) {*/
		fftw_real(Y, y);
		fft_real(x.data(), N);
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
		/*	} else {
		 auto a = fftw(Y);
		 auto b = FFT(X);
		 t1 += a;
		 t3 += a;
		 t2 += b;
		 t4 += b;
		 }*/
		for (int n = 0; n < N / 2 + 1; n++) {
			double x = X[n].real() - Y[n].real();
			double y = X[n].imag() - Y[n].imag();
			double err = sqrt(x * x + y * y);
			avg_err += err;
	//		printf("%e %e | %e %e | %e\n", X[n].real(), X[n].imag(), Y[n].real(), Y[n].imag(), err);
//			}
		}
		avg_err /= (21 * N);
		score *= cnt;
		score += t1 / (t2 + 1e-20);
		cnt++;
		score /= cnt;
		auto pfac = prime_factorization(N);
		std::string f;
		for (auto i = pfac.begin(); i != pfac.end(); i++) {
			f += "(" + std::to_string(i->first) + "^" + std::to_string(i->second) + ")";
		}
		printf("%i: %32s ", N, f.c_str());
		fflush(stdout);
		printf("| %e %e %e %e %e | %e\n", avg_err, t1, t2, t1 / (t2 + 1e-20), t4, score);
	}
	return 0;
}
