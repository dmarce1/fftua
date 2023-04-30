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
	feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

	constexpr int N = 64;
	std::vector<double> x(N);
	std::vector<complex<double>> y(N);
	for (int n = 0; n < N; n++) {
		y[n].imag() = 0.0;
		x[n] = y[n].real() = rand1();
	}
	fftw(y);
	fft_radix6step(x.data(), N);
	double err = abs(x[0] - y[0].real());
	std::vector<double> z(N);
	for (int n = 1; n < N - n; n++) {
		z[n] = y[n].real();
		z[N - n] = y[n].imag();
	}
	z[0] = y[0].real();
	z[N / 2] = y[N / 2].real();
	const int M = lround(sqrt(N));
	printf( "\n");
	for (int n = 0; n < M; n++) {
		printf("%4i | ", n);
		for (int m = 0; m < M; m++) {
			printf("%13.3e ", x[M * n + m]);
		}
		printf(" | ");
		for (int m = 0; m < M; m++) {
			printf("%13.3e ", z[M * n + m]);
		}
		printf(" | ");
		for (int m = 0; m < M; m++) {
			printf("%13.3e ", x[M * n + m] - z[M * n + m]);
		}
		printf(" \n");
		err += abs(x[n] - z[n]);

	}
	printf("%i | %e %e | %e %e\n", N / 2, x[N / 2], 0.0, y[N / 2].real(), y[N / 2].imag());
	err += abs(x[N / 2] - y[N / 2].real());
	printf("%e\n", err / N);
	//abort();
	timer tm3, tm4;
	double t3 = 0.0;
	double t4 = 0.0;
	std::vector<int> Ns;
	double score = 0.0;
	int cnt = 0;
	for (int N = 16; N <= 1024 * 1024 * 1024; N *= 4) {
		auto pfac = prime_factorization(N);
		{
			double avg_err = 0.0;
			double t1 = 0.0;
			double t2 = 0.0;
			for (int i = 0; i < 21; i++) {
				std::vector<double> x(N);
				std::vector<double> y(N);
				std::vector<complex<double>> X(N / 2 + 1);
				std::vector<complex<double>> Y(N / 2 + 1);
				for (int n = 0; n < N; n++) {
					x[n] = (y[n] = rand1());
				}
				if (i == 0) {
					fftw_real(Y, y);
					fft_radix6step(x.data(), N);
					//fht2fft(x.data(), N);
				} else {
					auto b = fftw_real(Y, y);
					timer tm;
					tm.start();
					fft_radix6step(x.data(), N);
					//	fht2fft(x.data(), N);
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
						//		printf("%16e %16e | %16e %16e | %16e %16e\n", X[n].real(), X[n].imag(), Y[n].real(), Y[n].imag(), X[n].real() - Y[n].real(), X[n].imag() - Y[n].imag());
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
			printf("R %c| %e %e %e %e %e | %e\n", (pfac.size() == 1 && pfac.begin()->second == 1) ? '*' : ' ', avg_err, t1, t2, t1 / (t2 + 1e-20), t4, score);
//				abort();
		}
		{
			double avg_err = 0.0;
			double t1 = 0.0;
			double t2 = 0.0;
			for (int i = 0; i < 21; i++) {
				std::vector<complex<double>> X(N);
				std::vector<complex<double>> Y(N);
				for (int n = 0; n < N; n++) {
					X[n].real() = (Y[n].real() = rand1());
					X[n].imag() = (Y[n].imag() = rand1());
				}
				if (i == 0) {
					fftw(Y);
					fft(X.data(), N);
				} else {
					auto b = fftw(Y);
					timer tm;
					tm.start();
					fft(X.data(), N);
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
						//		printf("%16e %16e | %16e %16e | %16e %16e\n", X[n].real(), X[n].imag(), Y[n].real(), Y[n].imag(), X[n].real() - Y[n].real(), X[n].imag() - Y[n].imag());
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
		}
		printf("\n");
	}
	return 0;
}
