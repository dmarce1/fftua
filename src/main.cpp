#include "fft.hpp"
#include "util.hpp"
#include "timer.hpp"
#include <cmath>

bool allowed_N(int N) {
	int k = 2;
	while (k <= SFFT_NMAX) {
		if (N % k == 0) {
			break;
		}
		k++;
	}
	if (k > SFFT_NMAX) {
		return false;
	}
	while (N > k) {
		N /= k;
		if (N % k != 0) {
			return false;
		}
	}
	return true;
}

int main(int argc, char **argv) {
	double t3 = 0.0;
	double t4 = 0.0;
	for (int N = 2; N <= 1024 * 1024; N++) {
		const auto pfac = prime_factorization(N);
		if (!allowed_N(N)) {
			continue;
		}
		double avg_err = 0.0;
		double t1 = 0.0;
		double t2 = 0.0;
		for (int i = 0; i < 25; i++) {
			std::vector<complex<double>> X(N);
			std::vector<complex<double>> Y(N);
			for (int n = 0; n < N; n++) {
				Y[n].real() = (X[n].real() = rand1());
				Y[n].imag() = (X[n].imag() = rand1());
			}
			if (i == 0) {
				fftw(Y);
				fft(X.data(), N);
			} else {
				timer tm;
				t1 += fftw(Y);
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
				//				printf("%e %e | %e %e | %e\n", X[n].real(), X[n].imag(), Y[n].real(), Y[n].imag(), err);
			}
		}
		avg_err /= (255 * N);
		std::string f;
		for (auto i = pfac.begin(); i != pfac.end(); i++) {
			f += "(" + std::to_string(i->first) + "^" + std::to_string(i->second) + ")";
		}
		t3 += t1;
		t4 += t2;
		printf("%i: %32s | %e | %e %e %e | %e %e %e\n", N, f.c_str(), avg_err, t1, t2, t1 / t2, t3, t4, t3 / t4);
	}
	return 0;
}

