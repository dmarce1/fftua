#include "fftu.hpp"
#include "util.hpp"

#include <cmath>

void fht_radix2(double* x, int N, bool scramble) {
	if (N <= 4) {
		sfht(x, N);
		return;
	}
	const int No4 = N >> 2;
	const auto& C = cos_twiddles(N);
	const auto& S = sin_twiddles(N);

	if (scramble) {
		int j = 1;
		for (int i = 1; i < N; i++) {
			if (i >= j) {
				std::swap(x[i - 1], x[j - 1]);
			}
			int k = No4;
			while (3 * k < j) {
				j -= 3 * k;
				k /= 4;
			}
			j += k;
		}
	}
	fht_radix2(x, N / 4, false);
	fht_radix2(x + N / 4, N / 4, false);
	fht_radix2(x + N / 2, N / 4, false);
	fht_radix2(x + 3 * N / 4, N / 4, false);
	const int L11 = 0;
	const int L12 = L11 + No4;
	const int L13 = L12 + No4;
	const int L14 = L13 + No4;
	const int L21 = No4 / 2;
	const int L22 = L21 + No4;
	const int L23 = L22 + No4;
	const int L24 = L23 + No4;
	auto T1 = x[L11] + x[L12];
	auto T2 = x[L11] - x[L12];
	auto T3 = x[L13] + x[L14];
	auto T4 = x[L13] - x[L14];
	x[L11] = T1 + T3;
	x[L12] = T1 - T3;
	x[L13] = T2 + T4;
	x[L14] = T2 - T4;
	T1 = x[L21];
	T2 = x[L22] * sqrt(2.0);
	T3 = x[L23];
	T4 = x[L24] * sqrt(2.0);
	x[L21] = T1 + T2 + T3;
	x[L22] = T1 - T3 + T4;
	x[L23] = T1 - T2 + T3;
	x[L24] = T1 - T3 - T4;
	for (int k = 1; k < N / 8; k++) {
		const auto C1 = C[k];
		const auto S1 = S[k];
		const auto C2 = C[2 * k];
		const auto S2 = S[2 * k];
		const auto C3 = C[3 * k];
		const auto S3 = S[3 * k];
		const int L11 = k;
		const int L12 = L11 + No4;
		const int L13 = L12 + No4;
		const int L14 = L13 + No4;
		const int L21 = No4 - k;
		const int L22 = L21 + No4;
		const int L23 = L22 + No4;
		const int L24 = L23 + No4;
		const auto T12 = x[L12] * C1 + x[L22] * S1;
		const auto T13 = x[L13] * C2 + x[L23] * S2;
		const auto T14 = x[L14] * C3 + x[L24] * S3;
		const auto T22 = x[L12] * S1 - x[L22] * C1;
		const auto T23 = x[L13] * S2 - x[L23] * C2;
		const auto T24 = x[L14] * S3 - x[L24] * C3;
		auto T1 = x[L21] + T23;
		auto T2 = x[L21] - T23;
		auto T3 = T22 + T24;
		auto T4 = T12 - T14;
		x[L24] = T2 - T3;
		x[L23] = T1 - T4;
		x[L22] = T2 + T3;
		x[L21] = T1 + T4;
		T1 = x[L11] + T13;
		T2 = x[L11] - T13;
		T3 = T24 - T22;
		T4 = T12 + T14;
		x[L14] = T2 - T3;
		x[L13] = T1 - T4;
		x[L12] = T2 + T3;
		x[L11] = T1 + T4;
	}
}
