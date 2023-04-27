#include "fftu.hpp"
#include "util.hpp"

#include <cmath>

void fht_radix22(double* x, int N) {
	if (N == 2) {
		const auto z0 = x[0];
		const auto z1 = x[1];
		x[0] = z0 + z1;
		x[1] = z0 - z1;
		return;
	} else if (N == 4) {
		const auto z0 = x[0];
		const auto z1 = x[1];
		const auto z2 = x[2];
		const auto z3 = x[3];
		const auto s0 = z0 + z1;
		const auto d0 = z0 - z1;
		const auto s1 = z2 + z3;
		const auto d1 = z2 - z3;
		x[0] = s0 + s1;
		x[1] = s0 - s1;
		x[2] = d0 + d1;
		x[3] = d0 - d1;
		return;
	}
	const auto& W = twiddles(N);
	std::vector<double> y(N);
	for (int n = 0; n < N / 2; n++) {
		y[n] = x[2 * n];
		y[n + N / 2] = x[2 * n + 1];
	}
	fht_radix2(y.data(), N / 2);
	fht_radix2(y.data() + N / 2, N / 2);
	for (int k = 0; k < N / 2; k++) {
		const auto& C = W[k].real();
		const auto& S = W[(N - k) % N].imag();
		const auto t = y[k + N / 2] * C + y[(N - k) % N] * S;
		x[k] = y[k] + t;
		x[k + N / 2] = y[k] - t;
	}
}

void fht_radix2(double* x, int N) {
	if (N == 2) {
		const auto z0 = x[0];
		const auto z1 = x[1];
		x[0] = z0 + z1;
		x[1] = z0 - z1;
		return;
	} else if (N == 4) {
		const auto z0 = x[0];
		const auto z1 = x[1];
		const auto z2 = x[2];
		const auto z3 = x[3];
		const auto s0 = z0 + z1;
		const auto d0 = z0 - z1;
		const auto s1 = z2 + z3;
		const auto d1 = z2 - z3;
		x[0] = s0 + s1;
		x[1] = s0 - s1;
		x[2] = d0 + d1;
		x[3] = d0 - d1;
		return;
	}
	const auto& W = twiddles(N);
	std::vector<double> y(N);
	for (int n = 0; n < N / 4; n++) {
		y[n] = x[4 * n];
		y[n + N / 4] = x[4 * n + 1];
		y[n + N / 2] = x[4 * n + 2];
		y[n + 3 * N / 4] = x[4 * n + 3];
	}
	fht_radix2(y.data(), N / 4);
	fht_radix2(y.data() + N / 4, N / 4);
	fht_radix2(y.data() + N / 2, N / 4);
	fht_radix2(y.data() + 3 * N / 4, N / 4);
	for (int k = 0; k < N / 4; k++) {

		const auto& C1 = W[k].real();
		const auto& S1 = W[mod(-k, N)].imag();
		const auto& C2 = W[2 * k].real();
		const auto& S2 = W[mod(-2 * k, N)].imag();
		const auto& C3 = W[3 * k].real();
		const auto& S3 = W[mod(-3 * k, N)].imag();
		const int kp0 = k;
		const int kp1 = k + N / 4;
		const int kp2 = k + N / 2;
		const int kp3 = k + 3 * N / 4;
		const int km = (N / 4 - k) % (N / 4);
		const int km0 = km;
		const int km1 = N / 4 + km;
		const int km2 = N / 2 + km;
		const int km3 = 3 * N / 4 + km;

		x[kp0] = y[k];
		x[kp1] = y[k];
		x[kp2] = y[k];
		x[kp3] = y[k];

		x[kp0] += +C1 * y[kp1] + S1 * y[km1];
		x[kp1] += -S1 * y[kp1] + C1 * y[km1];
		x[kp2] += -C1 * y[kp1] - S1 * y[km1];
		x[kp3] += +S1 * y[kp1] - C1 * y[km1];

		x[kp0] += +C2 * y[kp2] + S2 * y[km2];
		x[kp1] += -C2 * y[kp2] - S2 * y[km2];
		x[kp2] += +C2 * y[kp2] + S2 * y[km2];
		x[kp3] += -C2 * y[kp2] - S2 * y[km2];

		x[kp0] += +C3 * y[kp3] + S3 * y[km3];
		x[kp1] += +S3 * y[kp3] - C3 * y[km3];
		x[kp2] += -C3 * y[kp3] - S3 * y[km3];
		x[kp3] += -S3 * y[kp3] + C3 * y[km3];
	}
}
