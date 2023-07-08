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
void fft_2pow(double* x, int N);
void fft_2pow_complex(double* x, int N);
void fft_4step(double* X, int N);
void fft_inplace(double* x, int N);
void fft_inplace_real(double* x, int N);

extern "C" {
void fft_scramble(double* X, int N);
}
void test_twiddles();


extern "C" {
void twiddle_gen_next(__m256d* C, __m256d* S, void *ptr);
void twiddle_gen_init(void* ptr, int N);
void fft_simd_scramble(double*, int N);
void fft_transpose_hilo(double*, int, int);
void fft_recursive(double* X, const double* C, int N);
void dit_nr_recur(double* X,int N);
}


void test_twiddles() {
	void* gen;
	const int N = 1024*1024;
	posix_memalign(&gen, 32, 4 * 32 * ilogb(N) + 32);
	twiddle_gen_init(gen, N);
	__m256d C, S;
	for( int n = 0; n < 100; n++) {
	//	printf( "%e\n", *((double*)(gen)+n));
	}
//	return;
	for( int n = 1; n < N; n++) {
		twiddle_gen_next(&C, &S, gen);
		auto phi = 2.0 * M_PI * n / N;
//		auto phi = 2.0 * M_PI / (1 << (ilogb(N)-n));
		printf( "%i %e %e %e %e\n", n, C[0], S[0], cos(phi)-C[0], sin(-phi)-S[0]);
	}
	free(gen);
}

void test_time(double* x, int N) {
//	fft_simd_scramble(x, N);
//	fft_simd_scramble(x, N);
//	fft_scramble(x, N);
}

extern "C" {
int fft_bit_reverse(int, int);

}

int main(int argc, char **argv) {
	//test_twiddles();
	//return 0;
	constexpr int N = 256;
	timer tm;
	std::vector<double> x(N);
	for( int n = 0; n < N; n++) {
		x[n] = n;
	}
	tm.start();
//	fft_scramble(x.data(), N);
	tm.stop();
	//printf( "%e\n", tm.read());
	for( int n = 0; n < N; n++) {
		int i = n;
		int k = 0;
		for( int j = 0; j < ilogb(N); j++) {
			k <<= 1;
			k |= 1 & i;
			i >>= 1;
		}
	//	printf( "          .byte          %i %i\n", k, fft_bit_reverse(n,8));
	}
//	return 0;
//	feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
	timer tm3, tm4;
	double t3 = 0.0;
	double t4 = 0.0;
	std::vector<int> Ns;
	double score = 0.0;
	int cnt = 0;
	for (int N = 256; N <= 64*1024*1024; N *= 4) {
		auto pfac = prime_factorization(N);
		{
			double avg_err = 0.0;
			double t1 = 0.0;
			double t2 = 0.0;
			for (int i = 0; i < 1001; i++) {
				std::vector<double> x(N);
				std::vector<double> y(N);
				std::vector<complex<double>> X(N / 2 + 1);
				std::vector<complex<double>> Y(N / 2 + 1);
				for (int n = 0; n < N; n++) {
					x[n] = (y[n] = rand1());
				}
				x[1] = y[1] = 1.0;
//				x[0] = y[0] = 1.0;
				const auto& c = cos_twiddles(N);
				const auto& s = sin_twiddles(N);
				if (i == 0) {
					fftw_real(Y, y);
					dit_nr_recur(x.data(),N);
				} else {

					auto b = fftw_real(Y, y);
					timer tm;
					tm.start();
					dit_nr_recur(x.data(),  N);
					//test_time(x.data(), N);
					tm.stop();
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
					t2 += tm.read();
					t4 += tm.read();
					t1 += b;
					t3 += b;
					for (int n = 0; n < N / 2 + 1; n++) {
						double x = X[n].real() - Y[n].real();
						double y = X[n].imag() - Y[n].imag();
						double err = sqrt(x * x + y * y);
						avg_err += err;
			//			printf("%i: %16e %16e | %16e %16e | %16e %16e\n", n, X[n].real(), X[n].imag(), Y[n].real(), Y[n].imag(), X[n].real() - Y[n].real(), X[n].imag() - Y[n].imag());
					}
				}
			}
		//	abort();
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
			printf("R %c| %e %e %e %e %e | %e\n", (pfac.size() == 1 && pfac.begin()->second == 1) ? '*' : ' ', avg_err, t1, t2, t1 / (t2 + 1e-20), t4, score);
		}
		{
			continue;
			double avg_err = 0.0;
			double t1 = 0.0;
			double t2 = 0.0;
			for (int i = 0; i < 51; i++) {
				std::vector<complex<double>> X(N);
				std::vector<complex<double>> Y(N);
				for (int n = 0; n < N; n++) {
					X[n].real() = (Y[n].real() = rand1());
					X[n].imag() = (Y[n].imag() = rand1());
				}
			//	X[0].real() = (Y[1].real() = 1);
			//	X[0].imag() = (Y[0].imag() = 0);
				if (i == 0) {
					fftw(Y);
					fft_inplace((double*)X.data(), N);
				} else {
					auto b = fftw(Y);
					timer tm;
					tm.start();
					fft_inplace((double*)X.data(), N);
					tm.stop();
					t2 += tm.read();
					t4 += tm.read();
					t1 += b;
					t3 += b;
					for (int n = 0; n < N; n++) {
						double x = X[n].real() - Y[n].real();
						double y = X[n].imag() - Y[n].imag();
						double err = sqrt(x * x + y * y);
						avg_err += err;
			//			printf("%16e %16e | %16e %16e | %16e %16e\n", X[n].real(), X[n].imag(), Y[n].real(), Y[n].imag(), X[n].real() - Y[n].real(), X[n].imag() - Y[n].imag());
					}
				}
			}
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
			printf("C %c| %e %e %e %e %e | %e\n", (pfac.size() == 1 && pfac.begin()->second == 1) ? '*' : ' ', avg_err, t1, t2, t1 / (t2 + 1e-20), t4, score);
		//	abort();
		}
	}
	return 0;
}
