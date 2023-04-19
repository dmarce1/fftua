#include "fftu.hpp"
#include "util.hpp"
#include <cstring>

int compute_padding(int N) {
	static std::unordered_map<int, int> cache;
	auto iter = cache.find(N);
	if (iter == cache.end()) {
		int M = 2 * N - 1;
		bool done;
		do {
			done = true;
			auto factors = prime_factorization(M);
			const int two_pow = factors.begin()->first == 2 ? factors.begin()->second : 0;
			for (auto i = factors.begin(); i != factors.end(); i++) {
				if ((i->second > two_pow) || (i->first != 2 && i->first != 3 && i->first != 5)) {
					done = false;
					break;
				}
			}
			if (!done) {
				M++;
			}
		} while (!done);
		cache[N] = M;
		iter = cache.find(N);
	}
	return iter->second;
	//return 1 << (ilogb(2 * N - 1) + 1);
}

void fft_raders_indices(int* I, int N) {
	std::vector<int> J(N);
	const auto& gq = raders_gq(N);
	for (int n = 0; n < N - 1; n++) {
		J[n] = I[gq[n]];
	}
	J[N - 1] = I[0];
	std::memcpy(I, J.data(), N * sizeof(int));
	fft_indices(I, N - 1);
}

void fft_raders_padded_indices(int* I, int N) {
}

void fft_raders(complex<fft_simd4>* X, int N) {
	static thread_local std::unordered_map<int, std::vector<complex<fft_simd4>>>cache;
	std::vector<complex<fft_simd4>>& Y = cache[N];
	Y.resize(N - 1);
	const auto& gq = raders_gq(N);
	const auto& ginvq = raders_ginvq(N);
	const auto& tw = raders_twiddle(N, N - 1);
	complex<fft_simd4> xo = X[N - 1];
	complex<fft_simd4> x0 = xo;
	fft(X, N - 1);
	x0 += X[0];
	for (int q = 0; q < N - 1; q++) {
		Y[q] = X[q] * tw[q];
	}
	fft_scramble_inv(Y.data(), N - 1);
	fft(Y.data(), N - 1);
	for (int p = 0; p < N - 1; p++) {
		Y[p] += xo;
	}
	X[0] = x0;
	for (int p = 0; p < N - 1; p++) {
		X[ginvq[p]] = Y[p];
	}
}

void fft_raders_padded(complex<fft_simd4>* X, int N) {
	static thread_local std::unordered_map<int, std::vector<complex<fft_simd4>>>cache;
	std::vector<complex<fft_simd4>>& Y = cache[N];
	const int M = compute_padding(N);
	Y.resize(M);
	const auto& gq = raders_gq(N);
	const auto& ginvq = raders_ginvq(N);
	const auto& tw = raders_twiddle(N, M);
	complex<fft_simd4> xo = X[0];
	complex<fft_simd4> x0 = xo;
	for( int q = 0; q < N - 1; q++) {
		Y[q] = X[gq[q]];
		x0 += X[gq[q]];
	}
	for( int q = N - 1; q < M; q++) {
		Y[q] = fft_simd4(0.0);
	}
	fft_scramble(Y.data(), M);
	fft(Y.data(), M);
	for (int q = 0; q < M; q++) {
		Y[q] *= tw[q];
	}
	fft_scramble_inv(Y.data(), M);
	fft(Y.data(), M);
	for (int p = 0; p < N - 1; p++) {
		Y[p] += xo;
	}
	X[0] = x0;
	for (int p = 0; p < N - 1; p++) {
		X[ginvq[p]] = Y[p];
	}
}

void fft_raders(complex<double>* X, int N, bool padded) {
	static thread_local std::unordered_map<int, std::vector<complex<double>>>cache;
	std::vector<complex<double>>& Y = cache[N];
	const int M = padded ? compute_padding(N - 1) : N - 1;
	const int ysize = round_up(M, SIMD_SIZE);
	Y.resize(ysize);
	const auto& gq = raders_gq(N);
	const auto& ginvq = raders_ginvq(N);
	const auto& tw = raders_twiddle(N, M);
	fft_simd4* Z = (fft_simd4*) Y.data();
	complex<double> xo = X[0];
	complex<double> x0 = xo;
	fft_simd4 z0;
	for( int n = 1; n < SIMD_SIZE; n += 2) {
		z0[n] = xo.real();
	}
	for( int n = 0; n < SIMD_SIZE; n += 2) {
		z0[n] = xo.imag();
	}
	for (int n = 0; n < N - 1; n++) {
		Y[n] = X[gq[n]];
	}
	for (int n = N - 1; n < M; n++) {
		Y[n] = 0.0;
	}
	fft(Y.data(), M);
	x0 += Y[0];
	for (int q = 0; q < M; q += 2) {
		__m256d& y = *((__m256d*) &Y[q]);
		__m256d t = *((__m256d*) &tw[q]);
		y = mul(y, t);
		y = _mm256_permute_pd(y, 0x5);
	}
	fft(Y.data(), M);
	for (int p = 0; p < 2 * ysize / SIMD_SIZE; p++) {
		Z[p] += z0;
	}
	X[0] = x0;
	for (int p = 0; p < N - 1; p++) {
		const int pp = ginvq[p];
		X[pp].imag() = Y[p].real();
		X[pp].real() = Y[p].imag();
	}
}

