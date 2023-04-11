#include "fftu.hpp"
#include "fft.hpp"
#include "util.hpp"
#include <cstring>

void fft_stockham(int R, complex<fft_simd4>* X, int N) {
	auto scratch = create_scratch(N + R);
	auto* x0 = X;
	auto* x1 = scratch.data();
	auto* z = x1 + N;
	int NoR = N / R;
	int s = NoR;
	int sR = N;
	const auto& w = twiddles(N);
	while (s > 0) {
		for (int i = 0; i < s; i++) {
			int k = i;
			int j = i;
			for (int m = 0; m < NoR; m += s) {
				z[0] = x0[k];
				for (int l = 1; l < R; l++) {
					z[l] = x0[k + l * s] * w[m * l];
				}
		//		fft_complex_simd4((fft_simd4*) z, R);
				for (int l = 0; l < R; l++) {
					x1[j + l * NoR] = z[l];
				}
				k += sR;
				j += s;
			}
		}
		s /= R;
		sR /= R;
		auto* const tmp = x0;
		x0 = x1;
		x1 = tmp;
	}
	if (x0 == scratch.data()) {
		std::memcpy(x1, x0, N * sizeof(complex<fft_simd4> ));
	}
	destroy_scratch(std::move(scratch));
}
