#include "fftu.hpp"
#include "util.hpp"

void fft_odd(complex<double>* x, int N) {
	std::vector<complex<double>> y(N);
	const auto w = twiddles(N);
	for (int n = 0; n < N / 2; n++) {
		x[n] *= 2.0 * w[n];
	}
	fft(x, N);
}
