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
#include "fft.hpp"
#include "util.hpp"
#include "types.hpp"
#include "fftu.hpp"

void * operator new(std::size_t n) {
	void* memptr;
	posix_memalign(&memptr, 32, n);
	return memptr;
}

void operator delete(void * p) {
	free(p);
}

void *operator new[](std::size_t s) {
	void* memptr;
	posix_memalign(&memptr, 32, s);
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

int main(int argc, char **argv) {
	constexpr int w = 3;
	constexpr int N = 1 << w;
	timer tm3, tm4;
	double t3 = 0.0;
	double t4 = 0.0;
	for (int N = 8; N <= 1024 * 1024; N *= 2) {
		double avg_err = 0.0;
		double t1 = 0.0;
		double t2 = 0.0;
		for (int i = 0; i < 11; i++) {
			std::vector<complex<double>> X(N);
			std::vector<std::complex<double>> Y(N);
			for (int n = 0; n < N; n++) {
				Y[n].real(X[n].real() = rand1());
				Y[n].imag(X[n].imag() = rand1());
			}
			if (i == 0) {
				fftw(Y);
				FFT(X);
			} else {
				auto a = fftw(Y);
				auto b = FFT(X);
				t1 += a;
				t3 += a;
				t2 += b;
				t4 += b;
			}
			for (int n = 0; n < N; n++) {
				double x = X[n].real() - Y[n].real();
				double y = X[n].imag() - Y[n].imag();
				double err = sqrt(x * x + y * y);
				avg_err += err;
				//		printf("%e %e | %e %e | %e\n", X[n].real(), X[n].imag(), Y[n].real(), Y[n].imag(), err);
			}
		}
		avg_err /= (255 * N);
		auto pfac = prime_factorization(N);
		std::string f;
		for (auto i = pfac.begin(); i != pfac.end(); i++) {
			f += "(" + std::to_string(i->first) + "^" + std::to_string(i->second) + ")";
		}
		printf("%i: %32s | %e %e %e %e %e\n", N, f.c_str(), avg_err, t1, t2, t1 / t2, t4);
	}
	return 0;
}
