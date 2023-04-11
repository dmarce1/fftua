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

int swap_bits(int i, int b1, int b2) {
	int bit1 = (i >> b1) & 1;
	int bit2 = (i >> b2) & 1;
	i &= ~((1 << b1) | (1 << b2));
	i |= bit1 << b2;
	i |= bit2 << b1;
	return i;
}

int bit_reverse(int n, int w) {
	int i = n;
	int j = 0;
	for (int k = 0; k < w; k++) {
		j *= 2;
		j += i % 2;
		i /= 2;
	}
	return j;
}

void fft_self_sort(complex<double>* X, int N,  bool root = false) {
	if (N == 1) {
		return;
	}
	int nlev = ilogb(N);
	const auto& W = twiddles(N);
	for (int k = 0; k < N / 2; k++) {
		const auto z0 = X[k];
		const auto z1 = X[k + N / 2];
		X[k] = z0 + z1;
		X[k + N / 2] = z0 - z1;
		X[k + N / 2] *= W[k];
	}
	fft_self_sort(X, N / 2);
	fft_self_sort(X + N / 2,  N / 2);
	if (root) {
		for (int n = 0; n < N; n++) {
			const int l = bit_reverse(n, ilogb(N));
			if( n < l ) {
				std::swap(X[n], X[l]);
			}
		}
	}

}

int main(int argc, char **argv) {
	constexpr int w = 3;
	constexpr int N = 1 << w;
	timer tm3, tm4;
	for (int N = 8; N <= 64*1024 * 1024; N *= 2) {
		double avg_err = 0.0;
		double t1 = 0.0;
		double t2 = 0.0;
		for (int i = 0; i < 10; i++) {
			std::vector<complex<double>> X(N);
			std::vector<std::complex<double>> Y(N);
			for (int n = 0; n < N; n++) {
				Y[n].real(X[n].real() = rand1());
				Y[n].imag(X[n].imag() = rand1());
			}
			if (i == 0) {
				fftw(Y);
				fft(X.data(), N);
			} else {
				t1 += fftw(Y);
				timer tm;
				tm.start();
				fft(X.data(), N);
				tm.stop();
				t2 += tm.read();
			}
			for (int n = 0; n < N; n++) {
				double x = X[n].real() - Y[n].real();
				double y = X[n].imag() - Y[n].imag();
				double err = sqrt(x * x + y * y);
				avg_err += err;
				//	printf("%e %e | %e %e | %e\n", X[n].real(), X[n].imag(), Y[n].real(), Y[n].imag(), err);
			}
		}
		avg_err /= (255 * N);
		auto pfac = prime_factorization(N);
		std::string f;
		for (auto i = pfac.begin(); i != pfac.end(); i++) {
			f += "(" + std::to_string(i->first) + "^" + std::to_string(i->second) + ")";
		}
		printf("%i: %32s | %e %e %e %e\n", N, f.c_str(), avg_err, t1, t2, t1 / t2);
	}
	return 0;
}
