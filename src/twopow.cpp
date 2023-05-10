#include "fftu.hpp"
#include "util.hpp"
#include <memory>

int reverse_bits(int i, int N) {
	const auto w = ilogb(N);
	int j = 0;
	for (int k = 0; k < w; k++) {
		j <<= 1;
		j |= i & 1;
		i >>= 1;
	}
	return j;
}

int reverse_bits_conj(int i, int N) {
	const auto w = ilogb(N);
	int l = 0;
	int j = i;
	while (l < w) {
		int m = j % 4;
		if (m % 2 == 0) {
			l++;
			j >>= 1;
		} else {
			if (m % 4 == 3) {
				i -= 1 << (l + 2);
			}
			l += 2;
			j >>= 2;
		}
	}
	return reverse_bits(i, N);
}

template<class T>
void scramble(T* x, int N, int s = 1) {
	int N1m = N - 1;
	int j = 0;
	for (int i = 0; i < N1m; i++) {
		if (i > j) {
			std::swap(x[s * i], x[s * j]);
		}
		int k = N >> 1;
		while (k <= j) {
			j -= k;
			k >>= 1;
		}
		j += k;
	}
}

template<class T>
void scramble_complex(T* x, int N) {
	int N1m = N - 1;
	int j = 0;
	for (int i = 0; i < N1m; i++) {
		if (i > j) {
			std::swap(x[2 * i], x[2 * j]);
			std::swap(x[2 * i + 1], x[2 * j + 1]);
		}
		int k = N >> 1;
		while (k <= j) {
			j -= k;
			k >>= 1;
		}
		j += k;
	}
}

template<class T>
void scramble_complex4(T* x, int N) {
	int N1m = N - 1;
	int j = 0;
	for (int i = 0; i < N1m; i++) {
		if (i > j) {
			std::swap(x[2 * i], x[2 * j]);
			std::swap(x[2 * i + 1], x[2 * j + 1]);
		}
		int k = N / 4;
		while (3 * k <= j) {
			j -= 3 * k;
			k >>= 2;
		}
		j += k;
	}
}

template<class T>
void scramble_into(T* x, T* y, int N) {
	int N1m = N - 1;
	int j = 0;
	for (int i = 0; i < N1m; i++) {
		y[j] = x[i];
		int k = N >> 1;
		while (k <= j) {
			j -= k;
			k >>= 1;
		}
		j += k;
	}
	y[N1m] = x[N1m];
}

#include <numeric>

int int_pow(int a, int b) {
	int rc = 1;
	int apow = a;
	while (b) {
		if (b & 1) {
			rc *= apow;
		}
		b >>= 1;
		apow *= apow;
	}
	return rc;
}

void fft_4step1(fft_simd4** __restrict__ xre, fft_simd4** __restrict__ xim, int N, int M) {
	constexpr int N1 = 4;
	const int N2 = N / N1;
	const auto& W = twiddles(N);
	fft_simd4 C1, C2, C3, S1, S2, S3;
	int i0, i1, i2, i3;
	int j1, j2, j3;
	if (N == 1) {
		return;
	} else if (N == 2) {
		i0 = 0;
		i1 = i0 + 1;
		auto* U0R = xre[i0];
		auto* U0I = xim[i0];
		auto* U1R = xre[i1];
		auto* U1I = xim[i1];
		for (int m = 0; m < M; m++) {
			auto T0R = U0R[m];
			auto T0I = U0I[m];
			U0R[m] = T0R + U1R[m];
			U0I[m] = T0I + U1I[m];
			U1R[m] = T0R - U1R[m];
			U1I[m] = T0I - U1I[m];
		}
		return;
	}
	for (int n1 = 0; n1 < N1; n1++) {
		fft_4step1(xre + n1 * N2, xim + n1 * N2, N2, M);
	}
	const int chunk_size = std::min(M, 16);
	for (int chunk = 0; chunk < M / chunk_size; chunk++) {
		for (int k2 = 0; k2 < N2; k2++) {
			i0 = k2;
			i1 = i0 + N2;
			i2 = i1 + N2;
			i3 = i2 + N2;
			j1 = k2;
			j2 = 2 * j1;
			j3 = 3 * j1;
			C1 = W[j1].real();
			C2 = W[j2].real();
			C3 = W[j3].real();
			S1 = W[j1].imag();
			S2 = W[j2].imag();
			S3 = W[j3].imag();
			fft_simd4* __restrict__ U0R = xre[i0] + chunk * chunk_size;
			fft_simd4* __restrict__ U1R = xre[i2] + chunk * chunk_size;
			fft_simd4* __restrict__ U2R = xre[i1] + chunk * chunk_size;
			fft_simd4* __restrict__ U3R = xre[i3] + chunk * chunk_size;
			fft_simd4* __restrict__ U0I = xim[i0] + chunk * chunk_size;
			fft_simd4* __restrict__ U1I = xim[i2] + chunk * chunk_size;
			fft_simd4* __restrict__ U2I = xim[i1] + chunk * chunk_size;
			fft_simd4* __restrict__ U3I = xim[i3] + chunk * chunk_size;
			fft_simd4 T0R, T0I, T1R, T1I, T2R, T2I, T3R, T3I;
			for (int m = 0; m < chunk_size; m++) {
				T1R = U1R[m];
				U1R[m] = T1R * C1 - U1I[m] * S1;
				U1I[m] = T1R * S1 + U1I[m] * C1;
			}
			for (int m = 0; m < chunk_size; m++) {
				T2R = U2R[m];
				U2R[m] = T2R * C2 - U2I[m] * S2;
				U2I[m] = T2R * S2 + U2I[m] * C2;
			}
			for (int m = 0; m < chunk_size; m++) {
				auto T3R = U3R[m];
				U3R[m] = T3R * C3 - U3I[m] * S3;
				U3I[m] = T3R * S3 + U3I[m] * C3;
			}
			for (int m = 0; m < chunk_size; m++) {
				T0R = U0R[m];
				U0R[m] = T0R + U2R[m];
				U2R[m] = T0R - U2R[m];
			}
			for (int m = 0; m < chunk_size; m++) {
				T0I = U0I[m];
				U0I[m] = T0I + U2I[m];
				U2I[m] = T0I - U2I[m];
			}
			for (int m = 0; m < chunk_size; m++) {
				T1R = U1R[m];
				U1R[m] = T1R + U3R[m];
				U3R[m] = T1R - U3R[m];
			}
			for (int m = 0; m < chunk_size; m++) {
				T1I = U1I[m];
				U1I[m] = T1I + U3I[m];
				U3I[m] = T1I - U3I[m];
			}
			for (int m = 0; m < chunk_size; m++) {
				T0R = U0R[m];
				U0R[m] = T0R + U1R[m];
				U1R[m] = T0R - U1R[m];
			}
			for (int m = 0; m < chunk_size; m++) {
				T0I = U0I[m];
				U0I[m] = T0I + U1I[m];
				U1I[m] = T0I - U1I[m];
			}
			for (int m = 0; m < chunk_size; m++) {
				T2R = U2R[m];
				T2I = U2I[m];
				T3R = U3R[m];
				U2R[m] = T2R + U3I[m];
				U2I[m] = T2I - T3R;
				U3R[m] = T2R - U3I[m];
				U3I[m] = T2I + T3R;
			}
		}
	}
}

struct twiddle4 {
	double cr1;
	double cr2;
	double cr3;
	double ci1;
	double ci2;
	double ci3;
	double cr1n;
	double cr2n;
	double cr3n;
	double ci1n;
	double ci2n;
	double ci3n;
};

const std::vector<twiddle4>& twiddles4(int N) {
	using entry_type = std::shared_ptr<std::vector<twiddle4>>;
	static std::unordered_map<int, entry_type> cache;
	auto iter = cache.find(N);
	if (iter != cache.end()) {
		return *(iter->second);
	} else {
		std::vector<twiddle4> W(round_up(N, SIMD_SIZE));
		for (int n = 0; n < N; n++) {
			auto c1 = cos(-2.0 * M_PI * n / N);
			auto c2 = cos(-4.0 * M_PI * n / N);
			auto c3 = cos(-6.0 * M_PI * n / N);
			auto s1 = sin(-2.0 * M_PI * n / N);
			auto s2 = sin(-4.0 * M_PI * n / N);
			auto s3 = sin(-6.0 * M_PI * n / N);
			W[n].cr1 = c1;
			W[n].cr2 = c2;
			W[n].cr3 = c3 / c1;
			W[n].ci1 = s1 / c1;
			W[n].ci2 = s2 / c2;
			W[n].ci3 = s3 / c3;
			W[n].cr1n = -c1;
			W[n].cr2n = -c2;
			W[n].cr3n = -c3 / c1;
			W[n].ci1n = -s1 / c1;
			W[n].ci2n = -s2 / c2;
			W[n].ci3n = -s3 / c3;
		}
		cache[N] = std::make_shared<std::vector<twiddle4>>(std::move(W));
		return *(cache[N]);
	}
}

void fft_inplace(double* X, int N) {
	constexpr int N1 = 2;
	const auto& W = twiddles(N);
	int NHI = 1;
	int NLO = 1;
	int NMID = N / (N1 * N1);
	int N2 = N / N1;

	std::array<std::array<double, 2 * N1>, N1> U;
	std::array<double, N1> C;
	std::array<double, N1> S;

	const auto butterfly = [](std::array<double, 2 * N1>& U, const std::array<double, N1>& C, const std::array<double, N1>& S) {
		auto& u0r = U[0];
		auto& u0i = U[1];
		auto& u1r = U[2];
		auto& u1i = U[3];
		auto tmp = u1r;
		u1r = tmp * C[1] - u1i * S[1];
		u1i = tmp * S[1] + u1i * C[1];
		tmp = u0r;
		u0r += u1r;
		u1r = tmp - u1r;
		tmp = u0i;
		u0i += u1i;
		u1i = tmp - u1i;
	};

	while (NMID) {
	//	printf("NHI = %i NMID = %i NLO = %i\n", NHI, NMID, NLO);
		for (int nhi = 0; nhi < NHI; nhi++) {
			for (int nmid = 0; nmid < NMID; nmid++) {
				for (int nlo = 0; nlo < NLO; nlo++) {
					for (int na = 0; na < N1; na++) {
					//	printf("xi = ");
						for (int nb = 0; nb < N1; nb++) {
							const int i = nlo + NLO * (na + N1 * (nmid + NMID * (nb + N1 * nhi)));
							U[na][2 * nb] = X[2 * i];
							U[na][2 * nb + 1] = X[2 * i + 1];
			//				printf("%i ", i);
						}
					//	printf("\n");
					}
					//printf("ti = ");
					for (int na = 1; na < N1; na++) {
						const auto ti = na * nlo * NHI * N1 * NMID;
						const auto w = W[ti];
						C[na] = w.real();
						S[na] = w.imag();
					//	printf("%i ", ti);
					}
					//printf("\n");
					for (int na = 0; na < N1; na++) {
						butterfly(U[na], C, S);
					}
					for (int na = 0; na < N1; na++) {
						for (int nb = 0; nb < N1; nb++) {
							const int i = nlo + NLO * (na + N1 * (nmid + NMID * (nb + N1 * nhi)));
							X[2 * i] = U[nb][2 * na];
							X[2 * i + 1] = U[nb][2 * na + 1];
						}
					}
				}
			}
		}
		N2 /= N1;
		NHI *= N1;
		NLO *= N1;
		NMID /= N1 * N1;
	}
	if( (ilogb(N) / ilogb(N1)) % 2 == 0 ) {
		NHI /= N1;
	}
	while (NLO < N) {
		//printf("NHI = %i NLO = %i\n", NHI, NLO);
		for (int nhi = 0; nhi < NHI; nhi++) {
			for (int nlo = 0; nlo < NLO; nlo++) {
			//	printf("xi = ");
				for (int na = 0; na < N1; na++) {
					const int i = nlo + NLO * (na + N1 * nhi);
					U[0][2 * na] = X[2 * i];
					U[0][2 * na + 1] = X[2 * i + 1];
				//	printf("%i ", i);
				}
		//		printf("\n");
		//		printf("ti = ");
				for (int na = 1; na < N1; na++) {
					const int ti = na * nlo * NHI;
					const auto w = W[na * nlo * NHI];
					C[na] = w.real();
					S[na] = w.imag();
			//		printf("%i ", ti);
				}
		//		printf("\n");
				butterfly(U[0], C, S);
				for (int na = 0; na < N1; na++) {
					const int i = nlo + NLO * (na + N1 * nhi);
					X[2 * i] = U[0][2 * na];
					X[2 * i + 1] = U[0][2 * na + 1];
				}
			}
		}
		N2 *= N1;
		NHI /= N1;
		NLO *= N1;
	}
//	printf( "\n");
}

void fft_4step(double* X, int N) {
	const auto& W = twiddles(N);
	const int M = lround(sqrt(N));
	static std::vector<fft_simd4*> xre_ptrs;
	static std::vector<fft_simd4*> xim_ptrs;
	static std::vector<double> xre;
	static std::vector<double> xim;
	xim_ptrs.resize(N);
	xre_ptrs.resize(N);
	xre.resize(N);
	xim.resize(N);
	for (int n = 0; n < N; n++) {
		xre[n] = X[2 * n];
		xim[n] = X[2 * n + 1];
	}
	for (int m1 = 0; m1 < M; m1++) {
		scramble(xre.data() + m1, M, M);
		scramble(xim.data() + m1, M, M);
	}
	for (int m1 = 0; m1 < M; m1++) {
		xre_ptrs[m1] = ((fft_simd4*) (xre.data() + m1 * M));
		xim_ptrs[m1] = ((fft_simd4*) (xim.data() + m1 * M));
	}
	fft_4step1(xre_ptrs.data(), xim_ptrs.data(), M, M / SIMD_SIZE);
	for (int m1 = 0; m1 < M; m1++) {
		for (int m2 = 0; m2 < M; m2++) {
			const auto w = W[m1 * m2];
			auto& re = xre[m1 * M + m2];
			auto& im = xim[m1 * M + m2];
			auto tmp = re;
			re = tmp * w.real() - im * w.imag();
			im = tmp * w.imag() + im * w.real();
		}
	}
	for (int m1 = 0; m1 < M; m1++) {
		scramble(xre.data() + m1 * M, M);
		scramble(xim.data() + m1 * M, M);
	}
	for (int m1 = 0; m1 < M; m1++) {
		for (int m2 = 0; m2 < M; m2++) {
			const int i = m1 * M + m2;
			const int j = m2 * M + m1;
			if (i < j) {
				std::swap(xre[i], xre[j]);
				std::swap(xim[i], xim[j]);
			}
		}
	}
	fft_4step1(xre_ptrs.data(), xim_ptrs.data(), M, M / SIMD_SIZE);
	for (int n = 0; n < N; n++) {
		X[2 * n] = xre[n];
		X[2 * n + 1] = xim[n];
	}
}

void fft_width(double* X, int N) {
	constexpr int N1 = 4;

	const int npass = ilogb(N) / ilogb(N1);
	const auto& W = twiddles(N);

	const auto butterfly = [X, &W](int ki, int ti, int D, bool dif) {
		//	printf("pass = %i hi = %i lo = %i ki = %i ti = %i\n", pass, hi, lo, ki, ti);
			double T0R, T0I, T1R, T1I, T2R, T2I, T3R, T3I;
			double U0R, U0I, U1R, U1I, U2R, U2I, U3R, U3I;
			double C1, C2, C3, S1, S2, S3;
			int j1, j2, j3;
			int i0, i1, i2, i3, i4, i5, i6, i7;
			i0 = 2 * ki;
			i4 = 2 * ki + 1;
			i1 = i0 + 2 * D;
			i2 = i1 + 2 * D;
			i3 = i2 + 2 * D;
			i5 = i4 + 2 * D;
			i6 = i5 + 2 * D;
			i7 = i6 + 2 * D;
			if(!dif) {
				std::swap(i2, i1);
				std::swap(i5, i6);
			}
			U0R = X[i0];
			U1R = X[i1];
			U2R = X[i2];
			U3R = X[i3];
			U0I = X[i4];
			U1I = X[i5];
			U2I = X[i6];
			U3I = X[i7];
			if(!dif) {
				std::swap(i2, i1);
				std::swap(i5, i6);
			}
			j1 = ti;
			j2 = 2 * j1;
			j3 = 3 * j1;
			C1 = W[j1].real();
			C2 = W[j2].real();
			C3 = W[j3].real();
			S1 = W[j1].imag();
			S2 = W[j2].imag();
			S3 = W[j3].imag();
			if(!dif) {
				T1R = U1R;
				T2R = U2R;
				T3R = U3R;
				U1R = T1R * C1 - U1I * S1;
				U2R = T2R * C2 - U2I * S2;
				U3R = T3R * C3 - U3I * S3;
				U1I = T1R * S1 + U1I * C1;
				U2I = T2R * S2 + U2I * C2;
				U3I = T3R * S3 + U3I * C3;
			}
			T0R = U0R + U2R;
			T0I = U0I + U2I;
			T2R = U0R - U2R;

			T2I = U0I - U2I;
			T1R = U1R + U3R;

			T1I = U1I + U3I;
			T3R = U1R - U3R;
			T3I = U1I - U3I;

			U0R = T0R + T1R;
			U0I = T0I + T1I;
			U1R = T2R + T3I;
			U1I = T2I - T3R;
			U2R = T0R - T1R;
			U2I = T0I - T1I;

			U3R = T2R - T3I;
			U3I = T2I + T3R;

			if(dif) {
				T1R = U1R;
				T2R = U2R;
				T3R = U3R;
				U1R = T1R * C1 - U1I * S1;
				U2R = T2R * C2 - U2I * S2;
				U3R = T3R * C3 - U3I * S3;
				U1I = T1R * S1 + U1I * C1;
				U2I = T2R * S2 + U2I * C2;
				U3I = T3R * S3 + U3I * C3;
			}
			if(dif) {
				std::swap(i2, i1);
				std::swap(i5, i6);
			}
			X[i0] = U0R;
			X[i4] = U0I;
			X[i1] = U1R;
			X[i5] = U1I;
			X[i2] = U2R;
			X[i6] = U2I;
			X[i3] = U3R;
			X[i7] = U3I;
		};

	int nhi = 1;
	int nlo = N / N1;
	int np = 0;
	while (nlo > nhi / N1) {
		printf("nhi = %i nlo = %i\n", nhi, nlo);
		//		printf("pass = %i nhi = %i nlo = %i\n", pass, nhi, nlo);
		for (int hi = 0; hi < nhi; hi++) {
			for (int lo = 0; lo < nlo; lo++) {
				const int ki = hi * nlo * N1 + lo;
				const int ti = lo * nhi;
				butterfly(ki, ti, nlo, true);
			}
		}
		nhi *= N1;
		nlo /= N1;
		np++;
	}

	const auto M = lround(sqrt(N));
	for (int n1 = 0; n1 < M; n1++) {
		for (int k2 = 0; k2 < M; k2++) {
			const auto w = W[n1 * k2];
			auto& re = X[2 * (n1 * M + k2)];
			auto& im = X[2 * (n1 * M + k2) + 1];
			auto tmp = re;
			re = tmp * w.real() - im * w.imag();
			im = tmp * w.imag() + im * w.real();
		}
	}
	scramble_complex(X, N);
	nlo = 1;
	nhi = N / N1;
	while (nlo < N1 * nhi) {
		//	printf("nhi = %i nlo = %i\n", nhi, nlo);
//		printf("pass = %i nhi = %i nlo = %i\n", pass, nhi, nlo);
		for (int hi = 0; hi < nhi; hi++) {
			for (int lo = 0; lo < nlo; lo++) {
				const int ki = hi * nlo * N1 + lo;
				const int ti = lo * nhi;
				butterfly(ki, ti, nlo, false);
			}
		}
		nhi /= N1;
		nlo *= N1;
		np++;
	}
//	printf("%i\n", np);
}

#include <cstring>

void fft_width2(double* X0, int N) {
	constexpr int N1 = 4;
	constexpr int cache_size = 64;
	static std::vector<double> X1;
	X1.resize((N * 2));
	auto* X = X0;
	auto* Y = X1.data();
	int nepoch = ceil(log(N) / log(N1) / floor(log(cache_size) / log(N1)));
	const int nlev = ilogb(N) / ilogb(N1);
	while (nlev % nepoch != 0) {
		nepoch++;
	}
	const int npass = nlev / nepoch;
	int nbutter = std::min(cache_size, N);
	int ngroup = N / nbutter;
	while ((ilogb(ngroup) / ilogb(N1)) % (ilogb(nbutter) / ilogb(N1)) != 0) {
		nbutter /= N1;
		ngroup *= N1;
	}

//	printf("N1 = %i cache_size = %i nlev = %i nepoch = %i npass = %i ngroup = %i nbutter = %i\n", N1, cache_size, nlev, nepoch, npass, ngroup, nbutter);

	scramble_complex(X, N);
	const auto& W = twiddles(N);
	int nghi = ngroup;
	int nglo = 1;
	int digit = 0;
	for (int ei = 0; ei < nepoch; ei++) {
		for (int ghi = 0; ghi < nghi; ghi++) {
			for (int glo = 0; glo < nglo; glo++) {
				int this_digit = digit;
				int nbhi = nbutter / N1;
				int nblo = 1;
				for (int pi = 0; pi < npass; pi++) {
					const int tishft = int_pow(N1, nlev - this_digit - 1);
					for (int bhi = 0; bhi < nbhi; bhi++) {
						for (int blo = 0; blo < nblo; blo++) {
							const int k0 = blo + N1 * nblo * (bhi + nbhi * (glo + nglo * ghi));
							const int ti = (glo + nglo * blo) * tishft;
							//		printf("ei = %i ghi = %i glo = %i pi = %i bhi = %i blo = %i k0 = %i ti = %i nbhi = %i nblo = %i\n", ei, ghi, glo, pi, bhi, blo, k0, ti, nbhi, nblo);
							double T0R, T0I, T1R, T1I, T2R, T2I, T3R, T3I;
							double U0R, U0I, U1R, U1I, U2R, U2I, U3R, U3I;
							double C1, C2, C3, S1, S2, S3;
							int i0, i1, i2, i3, i4, i5, i6, i7;
							int j1, j2, j3;
							int D = nblo;
							i0 = 2 * k0;
							i4 = 2 * k0 + 1;
							i1 = i0 + 2 * D;
							i2 = i1 + 2 * D;
							i3 = i2 + 2 * D;
							i5 = i4 + 2 * D;
							i6 = i5 + 2 * D;
							i7 = i6 + 2 * D;
							U0R = X[i0];
							T1R = X[i2];
							T2R = X[i1];
							T3R = X[i3];
							U0I = X[i4];
							U1I = X[i6];
							U2I = X[i5];
							U3I = X[i7];
							j1 = ti;
							j2 = 2 * j1;
							j3 = 3 * j1;
							C1 = W[j1].real();
							C2 = W[j2].real();
							C3 = W[j3].real();
							S1 = W[j1].imag();
							S2 = W[j2].imag();
							S3 = W[j3].imag();
							U1R = T1R * C1 - U1I * S1;
							U2R = T2R * C2 - U2I * S2;
							U3R = T3R * C3 - U3I * S3;
							U1I = T1R * S1 + U1I * C1;
							U2I = T2R * S2 + U2I * C2;
							U3I = T3R * S3 + U3I * C3;
							T0R = U0R + U2R;
							T0I = U0I + U2I;
							T2R = U0R - U2R;
							T2I = U0I - U2I;
							T1R = U1R + U3R;
							T1I = U1I + U3I;
							T3R = U1R - U3R;
							T3I = U1I - U3I;
							U0R = T0R + T1R;
							U0I = T0I + T1I;
							U1R = T2R + T3I;
							U1I = T2I - T3R;
							U2R = T0R - T1R;
							U2I = T0I - T1I;
							U3R = T2R - T3I;
							U3I = T2I + T3R;
							X[i0] = U0R;
							X[i4] = U0I;
							X[i1] = U1R;
							X[i5] = U1I;
							X[i2] = U2R;
							X[i6] = U2I;
							X[i3] = U3R;
							X[i7] = U3I;
						}
					}
					nbhi /= N1;
					nblo *= N1;
					this_digit++;
					if (this_digit == nlev) {
						break;
					}
				}
			}
		}
		if (ei != nepoch - 1) {
			int nhi = nghi / nbutter;
			int nlo = nglo;
			for (int hi = 0; hi < nhi; hi++) {
				for (int n1 = 0; n1 < nbutter; n1++) {
					for (int n2 = 0; n2 < nbutter; n2++) {
						for (int lo = 0; lo < nlo; lo++) {
							const int i = n2 + nbutter * (lo + nlo * (n1 + nbutter * hi));
							const int j = n1 + nbutter * (lo + nlo * (n2 + nbutter * hi));
							Y[2 * i] = X[2 * j];
							Y[2 * i + 1] = X[2 * j + 1];
						}
					}
				}
			}
		} else {
			int N2 = nbutter;
			int N1 = nglo;
//			printf( "%i %i %i\n", N, N1, N2);
			for (int n2 = 0; n2 < N2; n2++) {
				for (int n1 = 0; n1 < N1; n1++) {
					const int i = n1 + N1 * n2;
					const int j = n2 + N2 * n1;
					//		printf( "%i %i %i %i %i %i %i\n", N, N1, N2, n1, n2, i, j);
					Y[2 * i] = X[2 * j];
					Y[2 * i + 1] = X[2 * j + 1];
				}
			}

		}
		std::swap(X, Y);
		nghi /= nbutter;
		nglo *= nbutter;
		digit += npass;
	}
	scramble_complex4(X, N);
	if (X != X0) {
		std::memcpy(X, Y, sizeof(double) * N);
	}
}

/*

 void fft_width(double* X, int N) {
 const auto& W = twiddles(N);
 constexpr int N1 = 4;
 double T0R, T0I, T1R, T1I, T2R, T2I, T3R, T3I;
 double U0R, U1R, U1I, U2R, U2I, U3R, U3I;
 double C1, C2, C3, S1, S2, S3;
 int i0, i1, i2, i3, i4, i5, i6, i7;
 int j1, j2, j3;
 int N2 = N / N1;
 int M1 = ilogb(N) - ilogb(N1);
 int M2 = M1 - ilogb(N1);

 std::vector<int> order(ilogb(N) / ilogb(N1));
 std::iota(order.begin(), order.end(), 0);
 const auto swap_indices = [N1, &order](int i, int m1, int m2) {
 if( i == 0 ) {
 std::swap(order[m1/ilogb(N1)], order[m2/ilogb(N1)]);
 for( int i = 0; i < order.size(); i++) {
 //			printf( "%i %i\n", i, order[i]);
 }
 //		printf( "\n");
 }
 const auto mask1 = (N1 - 1) << m1;
 const auto mask2 = (N1 - 1) << m2;
 auto index1 = (i & mask1) >> m1;
 auto index2 = (i & mask2) >> m2;
 i &= ~mask1;
 i &= ~mask2;
 index1 <<= m2;
 index2 <<= m1;
 return i | index1 | index2;
 };
 int iter = 0;
 int M = 1;
 int L = N / N1;
 while (L) {
 const auto& W = twiddles(M * N1);
 for (int l = 0; l < L; l++) {
 printf("%i %i %i\n", l, L, M);
 const int L0 = M * l;
 i0 = L0;
 i1 = i0 + N2;
 i2 = i1 + N2;
 i3 = i2 + N2;
 T0R = X[i0] + X[i2];
 T2R = X[i0] - X[i2];
 T1R = X[i1] + X[i3];
 T3R = X[i3] - X[i1];
 X[i0] = T0R + T1R;
 X[i1] = T2R;
 X[i2] = T0R - T1R;
 X[i3] = T3R;
 if (M * N1 >= 8) {
 i0 = L0 + N2 / 2;
 i1 = i0 + N2;
 i2 = i1 + N2;
 i3 = i2 + N2;
 printf("%i %i %i %i  : %i\n", i0, i1, i2, i3, N2);
 U0R = X[i0];
 U2R = X[i2];
 T1R = (X[i1] - X[i3]) * M_SQRT1_2;
 T3R = (X[i1] + X[i3]) * M_SQRT1_2;
 X[i0] = U0R + T1R;
 X[i1] = U0R - T1R;
 X[i2] = U2R - T3R;
 X[i3] = -U2R - T3R;
 }
 for (int k2 = 1; k2 < M / 2; k2++) {
 i0 = L0 + k2;
 i1 = i0 + N2;
 i2 = i1 + N2;
 i3 = i2 + N2;
 i4 = L0 + N2 - k2;
 i5 = i4 + N2;
 i6 = i5 + N2;
 i7 = i6 + N2;
 j1 = k2;
 j2 = 2 * j1;
 j3 = 3 * j1;
 printf("%i %i %i %i %i %i %i %i : %i\n", i0, i1, i2, i3, i4, i5, i6, i7, N2);
 C1 = W[j1].real();
 C2 = W[j2].real();
 C3 = W[j3].real();
 S1 = W[j1].imag();
 S2 = W[j2].imag();
 S3 = W[j3].imag();
 U1R = X[i1] * C1 - X[i5] * S1;
 U2R = X[i2] * C2 - X[i6] * S2;
 U3R = X[i3] * C3 - X[i7] * S3;
 U1I = X[i1] * S1 + X[i5] * C1;
 U2I = X[i2] * S2 + X[i6] * C2;
 U3I = X[i3] * S3 + X[i7] * C3;
 T0R = X[i0] + U2R;
 T0I = X[i4] + U2I;
 T2R = X[i0] - U2R;
 T2I = X[i4] - U2I;
 T1R = U1R + U3R;
 T1I = U1I + U3I;
 T3R = U3R - U1R;
 T3I = U1I - U3I;
 X[i0] = T0R + T1R;
 X[i1] = T2R + T3I;
 X[i2] = T1I - T0I;
 X[i3] = T3R - T2I;
 X[i4] = T2R - T3I;
 X[i5] = T0R - T1R;
 X[i6] = T2I + T3R;
 X[i7] = T0I + T1I;
 }
 }
 if (M2 >= 0) {
 for (int i = 0; i < N; i++) {
 int j = swap_indices(i, M1, M2);
 if (i < j) {
 std::swap(X[i], X[j]);
 }
 }
 }
 M2 -= ilogb(N1);
 M *= N1;
 L /= N1;
 }
 N2 = N / N1;
 for (int n1 = 0; n1 < N1; n1++) {
 int j = 0;
 for (int i = 0; i < N2 - 1; i++) {
 if (i > j) {
 //		std::swap(X[n1 * N2 + i], X[n1 * N2 + j]);
 }
 int k = N2 >> 1;
 while (k <= j) {
 j -= k;
 k >>= 1;
 }
 j += k;
 }
 }
 }*/

template<class T>
void fft_depth(T* X, int N) {
	constexpr int N1 = 4;
	if (N <= SFFT_NMAX) {
		scramble(X, N);
		sfft_real(X, N);
		return;
	}

	constexpr int N1o2 = N1 / 2;
	const int N2 = N / N1;
	const int No2 = N / 2;
	const int N2o2 = N2 / 2;

	static std::vector<int> D;
	if (!D.size()) {
		D.resize(N1);
		std::iota(D.begin(), D.end(), 0);
		scramble(D.data(), N1);
	}

	const auto& W = twiddles(N);

	for (int n1 = 0; n1 < N1; n1++) {
		fft_depth(&X[n1 * N2], N2);
	}

	std::array<complex<T>, N1> p;
	std::array<T, N1> q;
	for (int n1 = 0; n1 < N1; n1++) {
		q[n1] = X[N2 * D[n1]];
	}
	sfft_real(q.data(), N1);
	X[0] = q[0];
	X[No2] = q[N1o2];
	for (int k1 = 1; k1 < N1o2; k1++) {
		X[N2 * k1] = q[k1];
		X[N - N2 * k1] = q[N1 - k1];
	}
	if (N2o2) {
		for (int n1 = 0; n1 < N1; n1++) {
			q[n1] = X[N2 * D[n1] + N2o2];
		}
		sfft_skew(q.data(), N1);
		for (int k1 = 0; k1 < N1o2; k1++) {
			X[(N2 * k1 + N2 / 2)] = q[k1];
			X[N - (N2 * k1 + N2 / 2)] = q[N1 - k1 - 1];
		}
	}
	for (int k2 = 1; k2 < N2o2; k2++) {
		for (int n1 = 0; n1 < N1; n1++) {
			p[n1].real() = X[N2 * D[n1] + k2];
			p[n1].imag() = X[N2 * D[n1] - k2 + N2];
		}
		for (int n1 = 1; n1 < N1; n1++) {
			p[n1] *= W[n1 * k2];
		}
		sfft_complex((T*) p.data(), N1);
		for (int k1 = 0; k1 < N1o2; k1++) {
			const int k = N2 * k1 + k2;
			X[k] = p[k1].real();
			X[N - k] = p[k1].imag();
		}
		for (int k1 = N1o2; k1 < N1; k1++) {
			const int k = N2 * k1 + k2;
			X[N - k] = p[k1].real();
			X[k] = -p[k1].imag();
		}
	}
}

template<class T>
void fft_complex(T* X, T* Y, int N) {
	constexpr int N1 = 4;
	T T0R, T0I, T1R, T1I, T2R, T2I, T3R, T3I;
	T U0R, U0I, U1R, U1I, U2R, U2I, U3R, U3I;
	T C1, C2, C3, S1, S2, S3;
	int i0, i1, i2, i3;
	int j1, j2, j3;
	const int D = N / N1;
	int N2 = N / N1;
	int M1 = ilogb(N) - ilogb(N1);
	int M2 = M1 - ilogb(N1);
	const auto& W = twiddles(N);
	while (N2 >= 1) {
		const int M = N1 * N2;
		const int L = N / M;
		for (int k2 = 0; k2 < N / N1; k2++) {
			i0 = k2;
			i1 = i0 + D;
			i2 = i1 + D;
			i3 = i2 + D;
			U0R = X[i0];
			U1R = X[i1];
			U2R = X[i2];
			U3R = X[i3];
			U0I = Y[i0];
			U1I = Y[i1];
			U2I = Y[i2];
			U3I = Y[i3];
			T0R = U0R + U2R;
			T0I = U0I + U2I;
			T2R = U0R - U2R;
			T2I = U0I - U2I;
			T1R = U1R + U3R;
			T1I = U1I + U3I;
			T3R = U1R - U3R;
			T3I = U1I - U3I;
			U0R = T0R + T1R;
			U0I = T0I + T1I;
			U1R = T2R + T3I;
			U1I = T2I - T3R;
			U2R = T0R - T1R;
			U2I = T0I - T1I;
			U3R = T2R - T3I;
			U3I = T2I + T3R;
			j1 = L * (k2 % (M / N1));
			j2 = 2 * j1;
			j3 = 3 * j1;
			C1 = W[j1].real();
			C2 = W[j2].real();
			C3 = W[j3].real();
			S1 = W[j1].imag();
			S2 = W[j2].imag();
			S3 = W[j3].imag();
			T1R = U1R;
			T2R = U2R;
			T3R = U3R;
			U1R = T1R * C1 - U1I * S1;
			U2R = T2R * C2 - U2I * S2;
			U3R = T3R * C3 - U3I * S3;
			U1I = T1R * S1 + U1I * C1;
			U2I = T2R * S2 + U2I * C2;
			U3I = T3R * S3 + U3I * C3;
			X[i0] = U0R;
			Y[i0] = U0I;
			X[i1] = U1R;
			Y[i1] = U1I;
			X[i2] = U2R;
			Y[i2] = U2I;
			X[i3] = U3R;
			Y[i3] = U3I;
		}
		const int mask1 = (1 << M1) | (2 << M1);
		const int mask2 = (1 << M2) | (2 << M2);
		if (M2 >= 0) {
			for (int i = 0; i < N; i++) {
				int j = i;
				int p1 = (j & mask1) >> M1;
				int p2 = (j & mask2) >> M2;
				j &= ~mask1;
				j &= ~mask2;
				j |= p1 << M2;
				j |= p2 << M1;
				if (i < j) {
					std::swap(X[i], X[j]);
					std::swap(Y[i], Y[j]);
				}
			}
		}
		M2 -= ilogb(N1);
		N2 /= N1;
	}
	for (int n = 0; n < N1; n++) {
		int N2 = N / N1;
		auto* x = X + n * N / N1;
		auto* y = Y + n * N / N1;
		int N1m = N2 - 1;
		int j = 0;
		for (int i = 0; i < N1m; i++) {
			if (i > j) {
				std::swap(x[i], x[j]);
				std::swap(y[i], y[j]);
			}
			int k = N2 / 4;
			while (3 * k <= j) {
				j -= 3 * k;
				k >>= 2;
			}
			j += k;
		}
	}
}

template<class T>
void fft_complex2(T* X, T* Y, int N) {
	constexpr int N1 = 2;
	T T0R, T0I, T1R, T1I;
	T U0R, U0I, U1R, U1I;
	T C1, S1;
	int i0, i1;
	int j1;
	const int D = N / N1;
	int N2 = N / N1;
	int M1 = ilogb(N) - ilogb(N1);
	int M2 = M1 - ilogb(N1);
	while (N2 >= 1) {
		const int M = N1 * N2;
		for (int k2 = 0; k2 < N / N1; k2++) {
			const auto& W = twiddles(M);
			i0 = k2;
			i1 = i0 + D;
			U0R = X[i0];
			U1R = X[i1];
			U0I = Y[i0];
			U1I = Y[i1];
			T0R = U0R + U1R;
			T0I = U0I + U1I;
			T1R = U0R - U1R;
			T1I = U0I - U1I;
			U0R = T0R;
			U0I = T0I;
			U1R = T1R;
			U1I = T1I;
			j1 = k2 % (M / N1);
			C1 = W[j1].real();
			S1 = W[j1].imag();
			T1R = U1R;
			U1R = T1R * C1 - U1I * S1;
			U1I = T1R * S1 + U1I * C1;
			X[i0] = U0R;
			Y[i0] = U0I;
			X[i1] = U1R;
			Y[i1] = U1I;
		}
		const int mask1 = (1 << M1);
		const int mask2 = (1 << M2);
		if (M2 >= 0) {
			for (int i = 0; i < N; i++) {
				int j = i;
				int p1 = (j & mask1) >> M1;
				int p2 = (j & mask2) >> M2;
				j &= ~mask1;
				j &= ~mask2;
				j |= p1 << M2;
				j |= p2 << M1;
				if (i < j) {
					std::swap(X[i], X[j]);
					std::swap(Y[i], Y[j]);
				}
			}
		}
		M2 -= ilogb(N1);
		N2 /= N1;
	}
	for (int n = 0; n < 2; n++) {
		scramble(X + n * N / 2, N / 2);
		scramble(Y + n * N / 2, N / 2);
	}
}

template<class T>
void fft_dit(T* x, T* y, int N) {
	constexpr int N1 = 8;
	const int N2 = N / N1;
	T T0R, T0I, T1R, T1I, T2R, T2I, T3R, T3I;
	T T4R, T4I, T5R, T5I, T6R, T6I, T7R, T7I;
	T U0R, U0I, U1R, U1I, U2R, U2I, U3R, U3I;
	T U4R, U4I, U5R, U5I, U6R, U6I, U7R, U7I;
	T C1, C2, C3, C4, C5, C6, C7, S1, S2, S3, S4, S5, S6, S7;
	int i0, i1, i2, i3, i4, i5, i6, i7;
	int j1, j2, j3, j4, j5, j6, j7;

	if (N <= 1) {
		return;
	} else if (N == 2) {
		T0R = x[0];
		T0I = y[0];
		T1R = x[1];
		T1I = y[1];
		U0R = T0R + T1R;
		U0I = T0I + T1I;
		U1R = T0R - T1R;
		U1I = T0I - T1I;
		x[0] = U0R;
		x[1] = U1R;
		y[0] = U0I;
		y[1] = U1I;
		return;
	} else if (N == 4) {
		U0R = x[0];
		U1R = x[2];
		U2R = x[1];
		U3R = x[3];
		U0I = y[0];
		U1I = y[2];
		U2I = y[1];
		U3I = y[3];
		T0R = U0R + U2R;
		T0I = U0I + U2I;
		T2R = U0R - U2R;
		T2I = U0I - U2I;
		T1R = U1R + U3R;
		T1I = U1I + U3I;
		T3R = U1R - U3R;
		T3I = U1I - U3I;
		U0R = T0R + T1R;
		U0I = T0I + T1I;
		U1R = T2R + T3I;
		U1I = T2I - T3R;
		U2R = T0R - T1R;
		U2I = T0I - T1I;
		U3R = T2R - T3I;
		U3I = T2I + T3R;
		x[0] = U0R;
		y[0] = U0I;
		x[1] = U1R;
		y[1] = U1I;
		x[2] = U2R;
		y[2] = U2I;
		x[3] = U3R;
		y[3] = U3I;
	}

	for (int n1 = 0; n1 < N1; n1++) {
		fft_dit(x + n1 * N2, y + n1 * N2, N2);
	}
	const auto& W = twiddles(N);
	for (int k2 = 0; k2 < N2; k2++) {
		j1 = k2;
		j2 = 2 * k2;
		j3 = 3 * k2;
		j4 = 4 * k2;
		j5 = 5 * k2;
		j6 = 6 * k2;
		j7 = 7 * k2;
		C1 = W[j1].real();
		C2 = W[j2].real();
		C3 = W[j3].real();
		C4 = W[j4].real();
		C5 = W[j5].real();
		C6 = W[j6].real();
		C7 = W[j7].real();
		S1 = W[j1].imag();
		S2 = W[j2].imag();
		S3 = W[j3].imag();
		S4 = W[j4].imag();
		S5 = W[j5].imag();
		S6 = W[j6].imag();
		S7 = W[j7].imag();
		i0 = k2;
		i1 = i0 + N2;
		i2 = i1 + N2;
		i3 = i2 + N2;
		i4 = i3 + N2;
		i5 = i4 + N2;
		i6 = i5 + N2;
		i7 = i6 + N2;
		U0R = x[i0];
		U1R = x[i4];
		U2R = x[i2];
		U3R = x[i6];
		U4R = x[i1];
		U5R = x[i5];
		U6R = x[i3];
		U7R = x[i7];
		U0I = y[i0];
		U1I = y[i4];
		U2I = y[i2];
		U3I = y[i6];
		U4I = y[i1];
		U5I = y[i5];
		U6I = y[i3];
		U7I = y[i7];
		T1R = U1R;
		T2R = U2R;
		T3R = U3R;
		T4R = U4R;
		T5R = U5R;
		T6R = U6R;
		T7R = U7R;
		U1R = T1R * C1 - U1I * S1;
		U2R = T2R * C2 - U2I * S2;
		U3R = T3R * C3 - U3I * S3;
		U4R = T4R * C4 - U4I * S4;
		U5R = T5R * C5 - U5I * S5;
		U6R = T6R * C6 - U6I * S6;
		U7R = T7R * C7 - U7I * S7;
		U1I = T1R * S1 + U1I * C1;
		U2I = T2R * S2 + U2I * C2;
		U3I = T3R * S3 + U3I * C3;
		U4I = T4R * S4 + U4I * C4;
		U5I = T5R * S5 + U5I * C5;
		U6I = T6R * S6 + U6I * C6;
		U7I = T7R * S7 + U7I * C7;
		T0R = U0R + U4R;
		T0I = U0I + U4I;
		T2R = U0R - U4R;
		T2I = U0I - U4I;
		T1R = U2R + U6R;
		T1I = U2I + U6I;
		T3R = U2R - U6R;
		T3I = U2I - U6I;
		T4R = U1R + U5R;
		T4I = U1I + U5I;
		T6R = U1R - U5R;
		T6I = U1I - U5I;
		T5R = U3R + U7R;
		T5I = U3I + U7I;
		T7R = U3R - U7R;
		T7I = U3I - U7I;
		U0R = T0R + T1R;
		U0I = T0I + T1I;
		U2R = T2R + T3I;
		U2I = T2I - T3R;
		U4R = T0R - T1R;
		U4I = T0I - T1I;
		U6R = T2R - T3I;
		U6I = T2I + T3R;
		U1R = T4R + T5R;
		U1I = T4I + T5I;
		U3R = T6R + T7I;
		U3I = T6I - T7R;
		U5R = T4R - T5R;
		U5I = T4I - T5I;
		U7R = T6R - T7I;
		U7I = T6I + T7R;
		T3R = U3R;
		T5R = U5R;
		T7R = U7R;
		U3R = (T3R + U3I) * M_SQRT1_2;
		U3I = (-T3R + U3I) * M_SQRT1_2;
		U5R = U5I;
		U5I = -T5R;
		U7R = (-T7R + U7I) * M_SQRT1_2;
		U7I = (-T7R - U7I) * M_SQRT1_2;
		T0R = U0R;
		T0I = U0I;
		T2R = U2R;
		T2I = U2I;
		T4R = U4R;
		T4I = U4I;
		T6R = U6R;
		T6I = U6I;
		U0R = T0R + U1R;
		U0I = T0I + U1I;
		U1R = T0R - U1R;
		U1I = T0I - U1I;
		U2R = T2R + U3R;
		U2I = T2I + U3I;
		U3R = T2R - U3R;
		U3I = T2I - U3I;
		U4R = T4R + U5R;
		U4I = T4I + U5I;
		U5R = T4R - U5R;
		U5I = T4I - U5I;
		U6R = T6R + U7R;
		U6I = T6I + U7I;
		U7R = T6R - U7R;
		U7I = T6I - U7I;
		x[i0] = U0R;
		y[i0] = U0I;
		x[i1] = U2R;
		y[i1] = U2I;
		x[i2] = U4R;
		y[i2] = U4I;
		x[i3] = U6R;
		y[i3] = U6I;
		x[i4] = U1R;
		y[i4] = U1I;
		x[i5] = U3R;
		y[i5] = U3I;
		x[i6] = U5R;
		y[i6] = U5I;
		x[i7] = U7R;
		y[i7] = U7I;
	}
}

void fft_2pow_complex(double* x, int N) {
	std::vector<double> re(N), im(N);
	for (int n = 0; n < N; n++) {
		re[n] = x[2 * n];
		im[n] = x[2 * n + 1];
	}

	scramble(re.data(), N);
	scramble(im.data(), N);
	fft_dit(re.data(), im.data(), N);

	for (int n = 0; n < N; n++) {
		x[2 * n] = re[n];
		x[2 * n + 1] = im[n];
	}
}

template<class T>
void fft_dit_real(T* x, int N) {
	constexpr int N1 = 8;
	const int N2 = N / N1;
	T T0R, T0I, T1R, T1I, T2R, T2I, T3R, T3I;
	T T4R, T4I, T5R, T5I, T6R, T6I, T7R, T7I;
	T U0R, U0I, U1R, U1I, U2R, U2I, U3R, U3I;
	T U4R, U4I, U5R, U5I, U6R, U6I, U7R, U7I;
	T C1, C2, C3, C4, C5, C6, C7, S1, S2, S3, S4, S5, S6, S7;
	int i0, i1, i2, i3, i4, i5, i6, i7;
	int i8, i9, i10, i11, i12, i13, i14, i15;
	int j1, j2, j3, j4, j5, j6, j7;

	if (N <= 1) {
		return;
	} else if (N == 2) {
		T0R = x[0];
		T1R = x[1];
		U0R = T0R + T1R;
		U1R = T0R - T1R;
		x[0] = U0R;
		x[1] = U1R;
		return;
	} else if (N == 4) {
		U0R = x[0];
		U1R = x[2];
		U2R = x[1];
		U3R = x[3];
		T0R = U0R + U2R;
		T2R = U0R - U2R;
		T1R = U1R + U3R;
		T3R = U1R - U3R;
		U0R = T0R + T1R;
		U1R = T2R;
		U1I = -T3R;
		U2R = T0R - T1R;
		x[0] = U0R;
		x[1] = U1R;
		x[2] = U2R;
		x[3] = U1I;
		return;
	}

	for (int n1 = 0; n1 < N1; n1++) {
		fft_dit_real(x + n1 * N2, N2);
	}
	const auto& W = twiddles(N);
	i0 = 0;
	i1 = i0 + N2;
	i2 = i1 + N2;
	i3 = i2 + N2;
	i4 = i3 + N2;
	i5 = i4 + N2;
	i6 = i5 + N2;
	i7 = i6 + N2;
	U0R = x[i0];
	U1R = x[i4];
	U2R = x[i2];
	U3R = x[i6];
	U4R = x[i1];
	U5R = x[i5];
	U6R = x[i3];
	U7R = x[i7];
	T0R = U0R + U4R;
	T2R = U0R - U4R;
	T1R = U2R + U6R;
	T3R = U2R - U6R;
	T4R = U1R + U5R;
	T6R = U1R - U5R;
	T5R = U3R + U7R;
	T7R = U3R - U7R;
	U0R = T0R + T1R;
	U2R = T2R;
	U2I = -T3R;
	U4R = T0R - T1R;
	U1R = T4R + T5R;
	U3R = T6R;
	U3I = -T7R;
	U5R = T4R - T5R;
	T3R = U3R;
	T5R = U5R;
	U3R = (T3R + U3I) * M_SQRT1_2;
	U3I = (-T3R + U3I) * M_SQRT1_2;
	U5I = -T5R;
	T0R = U0R;
	T2R = U2R;
	T2I = U2I;
	T4R = U4R;
	T4I = U4I;
	T6R = U2R;
	T6I = -U2I;
	U0R = T0R + U1R;
	U1R = T0R - U1R;
	U2R = T2R + U3R;
	U2I = T2I + U3I;
	U4R = T4R;
	U4I = U5I;
	U6R = T6R - U3R;
	U6I = T6I + U3I;
	x[i0] = U0R;
	x[i1] = U2R;
	x[i2] = U4R;
	x[i3] = U6R;
	x[i4] = U1R;
	x[i5] = U6I;
	x[i6] = U4I;
	x[i7] = U2I;
	if (N2 / 2) {
		constexpr double c1 = cos((1.0 / 8.0) * M_PI);
		constexpr double c3 = cos((3.0 / 8.0) * M_PI);
		i0 = N2 / 2;
		i1 = i0 + N2;
		i2 = i1 + N2;
		i3 = i2 + N2;
		i4 = i3 + N2;
		i5 = i4 + N2;
		i6 = i5 + N2;
		i7 = i6 + N2;
		U0R = x[i0];
		U1R = x[i4];
		U2R = x[i2];
		U3R = x[i6];
		U4R = x[i1];
		U5R = x[i5];
		U6R = x[i3];
		U7R = x[i7];
		T1R = U1R;
		T2R = U2R;
		T3R = U3R;
		T5R = U5R;
		T6R = U6R;
		T7R = U7R;
		U1R = T1R * c1;
		U3R = T3R * c3;
		U5R = T5R * c3;
		U7R = T7R * c1;
		U1I = T1R * c3;
		U3I = T3R * c1;
		U5I = T5R * c1;
		U7I = T7R * c3;
		T0R = U0R;
		T0I = U4R;
		T1R = (U2R - U6R) * M_SQRT1_2;
		T3R = (U2R + U6R) * M_SQRT1_2;
		T4R = U1R - U5R;
		T4I = -U1I - U5I;
		T6R = U1R + U5R;
		T6I = -U1I + U5I;
		T5R = U3R - U7R;
		T5I = -U3I - U7I;
		T7R = U3R + U7R;
		T7I = -U3I + U7I;
		U0R = T0R + T1R;
		U0I = -T0I - T3R;
		U2R = T0R - T1R;
		U2I = T0I - T3R;
		U6R = T0R + T1R;
		U6I = T0I + T3R;
		U1R = T4R + T5R;
		U1I = T4I + T5I;
		U3R = T6R + T7I;
		U3I = T6I - T7R;
		U5R = T4R - T5R;
		U5I = T4I - T5I;
		U7R = T6R - T7I;
		U7I = T6I + T7R;
		T3R = U3R;
		T5R = U5R;
		T7R = U7R;
		U3R = (T3R + U3I) * M_SQRT1_2;
		U3I = (-T3R + U3I) * M_SQRT1_2;
		U5R = U5I;
		U5I = -T5R;
		U7R = (-T7R + U7I) * M_SQRT1_2;
		U7I = (-T7R - U7I) * M_SQRT1_2;
		U0R = U0R + U1R;
		U0I = U0I + U1I;
		U4R = U2R + U5R;
		U4I = -U2I + U5I;
		U2R = U2R + U3R;
		U2I = U2I + U3I;
		U6R = U6R + U7R;
		U6I = U6I + U7I;
		x[i0] = U0R;
		x[i7] = U0I;
		x[i1] = U2R;
		x[i6] = U2I;
		x[i2] = U4R;
		x[i5] = U4I;
		x[i3] = U6R;
		x[i4] = U6I;
	}
	for (int k2 = 1; k2 < N2 / 2; k2++) {
		j1 = k2;
		j2 = 2 * k2;
		j3 = 3 * k2;
		j4 = 4 * k2;
		j5 = 5 * k2;
		j6 = 6 * k2;
		j7 = 7 * k2;
		C1 = W[j1].real();
		C2 = W[j2].real();
		C3 = W[j3].real();
		C4 = W[j4].real();
		C5 = W[j5].real();
		C6 = W[j6].real();
		C7 = W[j7].real();
		S1 = W[j1].imag();
		S2 = W[j2].imag();
		S3 = W[j3].imag();
		S4 = W[j4].imag();
		S5 = W[j5].imag();
		S6 = W[j6].imag();
		S7 = W[j7].imag();
		i0 = k2;
		i1 = i0 + N2;
		i2 = i1 + N2;
		i3 = i2 + N2;
		i4 = i3 + N2;
		i5 = i4 + N2;
		i6 = i5 + N2;
		i7 = i6 + N2;
		i8 = N2 - k2;
		i9 = i8 + N2;
		i10 = i9 + N2;
		i11 = i10 + N2;
		i12 = i11 + N2;
		i13 = i12 + N2;
		i14 = i13 + N2;
		i15 = i14 + N2;
		U0R = x[i0];
		U0I = x[i8];
		U1R = x[i4];
		U1I = x[i12];
		U2R = x[i2];
		U2I = x[i10];
		U3R = x[i6];
		U3I = x[i14];
		U4R = x[i1];
		U4I = x[i9];
		U5R = x[i5];
		U5I = x[i13];
		U6R = x[i3];
		U6I = x[i11];
		U7R = x[i7];
		U7I = x[i15];
		T1R = U1R;
		T2R = U2R;
		T3R = U3R;
		T4R = U4R;
		T5R = U5R;
		T6R = U6R;
		T7R = U7R;
		U1R = T1R * C1 - U1I * S1;
		U3R = T3R * C3 - U3I * S3;
		U5R = T5R * C5 - U5I * S5;
		U7R = T7R * C7 - U7I * S7;
		U1I = T1R * S1 + U1I * C1;
		U3I = T3R * S3 + U3I * C3;
		U5I = T5R * S5 + U5I * C5;
		U7I = T7R * S7 + U7I * C7;
		U2R = T2R * C2 - U2I * S2;
		U4R = T4R * C4 - U4I * S4;
		U6R = T6R * C6 - U6I * S6;
		U2I = T2R * S2 + U2I * C2;
		U4I = T4R * S4 + U4I * C4;
		U6I = T6R * S6 + U6I * C6;
		T0R = U0R + U4R;
		T0I = U0I + U4I;
		T2R = U0R - U4R;
		T2I = U0I - U4I;
		T1R = U2R + U6R;
		T1I = U2I + U6I;
		T3R = U2R - U6R;
		T3I = U2I - U6I;
		T4R = U1R + U5R;
		T4I = U1I + U5I;
		T6R = U1R - U5R;
		T6I = U1I - U5I;
		T5R = U3R + U7R;
		T5I = U3I + U7I;
		T7R = U3R - U7R;
		T7I = U3I - U7I;
		U0R = T0R + T1R;
		U0I = T0I + T1I;
		U2R = T2R + T3I;
		U2I = T2I - T3R;
		U4R = T0R - T1R;
		U4I = T0I - T1I;
		U6R = T2R - T3I;
		U6I = T2I + T3R;
		U1R = T4R + T5R;
		U1I = T4I + T5I;
		U3R = T6R + T7I;
		U3I = T6I - T7R;
		U5R = T4R - T5R;
		U5I = T4I - T5I;
		U7R = T6R - T7I;
		U7I = T6I + T7R;
		T3R = U3R;
		T5R = U5R;
		T7R = U7R;
		U3R = (T3R + U3I) * M_SQRT1_2;
		U3I = (-T3R + U3I) * M_SQRT1_2;
		U5R = U5I;
		U5I = -T5R;
		U7R = (-T7R + U7I) * M_SQRT1_2;
		U7I = (-T7R - U7I) * M_SQRT1_2;
		T0R = U0R;
		T0I = U0I;
		T2R = U2R;
		T2I = U2I;
		T4R = U4R;
		T4I = U4I;
		T6R = U6R;
		T6I = U6I;
		U0R = T0R + U1R;
		U0I = T0I + U1I;
		U1R = T0R - U1R;
		U1I = T0I - U1I;
		U2R = T2R + U3R;
		U2I = T2I + U3I;
		U3R = T2R - U3R;
		U3I = T2I - U3I;
		U4R = T4R + U5R;
		U4I = T4I + U5I;
		U5R = T4R - U5R;
		U5I = T4I - U5I;
		U6R = T6R + U7R;
		U6I = T6I + U7I;
		U7R = T6R - U7R;
		U7I = T6I - U7I;
		x[i0] = U0R;
		x[i15] = U0I;
		x[i1] = U2R;
		x[i14] = U2I;
		x[i2] = U4R;
		x[i13] = U4I;
		x[i3] = U6R;
		x[i12] = U6I;
		x[i11] = U1R;
		x[i4] = -U1I;
		x[i10] = U3R;
		x[i5] = -U3I;
		x[i9] = U5R;
		x[i6] = -U5I;
		x[i8] = U7R;
		x[i7] = -U7I;
	}
}

template<class T>
void fft_dit_real_split1(T* x, int N) {
	constexpr int N1 = 8;
	const int N2 = N / N1;
	T T0R, T0I, T1R, T1I, T2R, T2I, T3R, T3I;
	T T4R, T4I, T5R, T5I, T6R, T6I, T7R, T7I;
	T U0R, U0I, U1R, U1I, U2R, U2I, U3R, U3I;
	T U4R, U4I, U5R, U5I, U6R, U6I, U7R, U7I;
	T C1, C2, C3, C5, C6, C7, S1, S2, S3, S5, S6, S7;
	int i0, i1, i2, i3, i4, i5, i6, i7;
	int i8, i9, i10, i11, i12, i13, i14, i15;
	int j1, j2, j3, j5, j6, j7;

	if (N <= 1) {
		return;
	} else if (N == 2) {
		T0R = x[0];
		T1R = x[1];
		U0R = T0R + T1R;
		U1R = T0R - T1R;
		x[0] = U0R;
		x[1] = U1R;
		return;
	} else if (N == 4) {
		U0R = x[0];
		U1R = x[2];
		U2R = x[1];
		U3R = x[3];
		T0R = U0R + U2R;
		T2R = U0R - U2R;
		T1R = U1R + U3R;
		T3R = U1R - U3R;
		U0R = T0R + T1R;
		U1R = T2R;
		U1I = -T3R;
		U2R = T0R - T1R;
		x[0] = U0R;
		x[1] = U1R;
		x[2] = U2R;
		x[3] = U1I;
		return;
	}
	fft_dit_real_split1(x, 2 * N2);
	for (int n1 = 2; n1 < N1; n1++) {
		fft_dit_real_split1(x + n1 * N2, N2);
	}
	const auto& W = twiddles(N);
	i0 = 0;
	i1 = i0 + N2;
	i2 = i1 + N2;
	i3 = i2 + N2;
	i4 = i3 + N2;
	i5 = i4 + N2;
	i6 = i5 + N2;
	i7 = i6 + N2;
	U0R = x[i0];
	U1R = x[i4];
	U2R = x[i2];
	U3R = x[i6];
	U4R = x[i1];
	U5R = x[i5];
	U6R = x[i3];
	U7R = x[i7];
	T0R = U0R;
	T2R = U4R;
	T1R = U2R + U6R;
	T3R = U2R - U6R;
	T4R = U1R + U5R;
	T6R = U1R - U5R;
	T5R = U3R + U7R;
	T7R = U3R - U7R;
	U0R = T0R + T1R;
	U2R = T2R;
	U2I = -T3R;
	U4R = T0R - T1R;
	U1R = T4R + T5R;
	U3R = T6R;
	U3I = -T7R;
	U5R = T4R - T5R;
	T3R = U3R;
	T5R = U5R;
	U3R = (T3R + U3I) * M_SQRT1_2;
	U3I = (-T3R + U3I) * M_SQRT1_2;
	U5I = -T5R;
	T0R = U0R;
	T2R = U2R;
	T2I = U2I;
	T4R = U4R;
	T4I = U4I;
	T6R = U2R;
	T6I = -U2I;
	U0R = T0R + U1R;
	U1R = T0R - U1R;
	U2R = T2R + U3R;
	U2I = T2I + U3I;
	U4R = T4R;
	U4I = U5I;
	U6R = T6R - U3R;
	U6I = T6I + U3I;
	x[i0] = U0R;
	x[i1] = U2R;
	x[i2] = U4R;
	x[i3] = U6R;
	x[i4] = U1R;
	x[i5] = U6I;
	x[i6] = U4I;
	x[i7] = U2I;
	if (N2 / 2) {
		constexpr double c1 = cos((1.0 / 8.0) * M_PI);
		constexpr double c3 = cos((3.0 / 8.0) * M_PI);
		i0 = N2 / 2;
		i1 = i0 + N2;
		i2 = i1 + N2;
		i3 = i2 + N2;
		i4 = i3 + N2;
		i5 = i4 + N2;
		i6 = i5 + N2;
		i7 = i6 + N2;
		U0R = x[i0];
		U1R = x[i4];
		U2R = x[i2];
		U3R = x[i6];
		U4R = x[i1];
		U5R = x[i5];
		U6R = x[i3];
		U7R = x[i7];
		T1R = U1R;
		T2R = U2R;
		T3R = U3R;
		T5R = U5R;
		T6R = U6R;
		T7R = U7R;
		U1R = T1R * c1;
		U3R = T3R * c3;
		U5R = T5R * c3;
		U7R = T7R * c1;
		U1I = T1R * c3;
		U3I = T3R * c1;
		U5I = T5R * c1;
		U7I = T7R * c3;
		T0R = U0R;
		T0I = -U4R;
		T1R = (U2R - U6R) * M_SQRT1_2;
		T3R = (U2R + U6R) * M_SQRT1_2;
		T4R = U1R - U5R;
		T4I = -U1I - U5I;
		T6R = U1R + U5R;
		T6I = -U1I + U5I;
		T5R = U3R - U7R;
		T5I = -U3I - U7I;
		T7R = U3R + U7R;
		T7I = -U3I + U7I;
		U0R = T0R + T1R;
		U0I = -T0I - T3R;
		U2R = T0R - T1R;
		U2I = T0I - T3R;
		U6R = T0R + T1R;
		U6I = T0I + T3R;
		U1R = T4R + T5R;
		U1I = T4I + T5I;
		U3R = T6R + T7I;
		U3I = T6I - T7R;
		U5R = T4R - T5R;
		U5I = T4I - T5I;
		U7R = T6R - T7I;
		U7I = T6I + T7R;
		T3R = U3R;
		T5R = U5R;
		T7R = U7R;
		U3R = (T3R + U3I) * M_SQRT1_2;
		U3I = (-T3R + U3I) * M_SQRT1_2;
		U5R = U5I;
		U5I = -T5R;
		U7R = (-T7R + U7I) * M_SQRT1_2;
		U7I = (-T7R - U7I) * M_SQRT1_2;
		U0R = U0R + U1R;
		U0I = U0I + U1I;
		U4R = U2R + U5R;
		U4I = -U2I + U5I;
		U2R = U2R + U3R;
		U2I = U2I + U3I;
		U6R = U6R + U7R;
		U6I = U6I + U7I;
		x[i0] = U0R;
		x[i7] = U0I;
		x[i1] = U2R;
		x[i6] = U2I;
		x[i2] = U4R;
		x[i5] = U4I;
		x[i3] = U6R;
		x[i4] = U6I;
	}
	for (int k2 = 1; k2 < N2 / 2; k2++) {
		j1 = k2;
		j2 = 2 * k2;
		j3 = 3 * k2;
		j5 = 5 * k2;
		j6 = 6 * k2;
		j7 = 7 * k2;
		C1 = W[j1].real();
		C2 = W[j2].real();
		C3 = W[j3].real();
		C5 = W[j5].real();
		C6 = W[j6].real();
		C7 = W[j7].real();
		S1 = W[j1].imag();
		S2 = W[j2].imag();
		S3 = W[j3].imag();
		S5 = W[j5].imag();
		S6 = W[j6].imag();
		S7 = W[j7].imag();
		i0 = k2;
		i1 = i0 + N2;
		i2 = i1 + N2;
		i3 = i2 + N2;
		i4 = i3 + N2;
		i5 = i4 + N2;
		i6 = i5 + N2;
		i7 = i6 + N2;
		i8 = N2 - k2;
		i9 = i8 + N2;
		i10 = i9 + N2;
		i11 = i10 + N2;
		i12 = i11 + N2;
		i13 = i12 + N2;
		i14 = i13 + N2;
		i15 = i14 + N2;
		U0R = x[i0];
		U4R = x[i8];
		U4I = x[i1];
		U0I = x[i9];
		U1R = x[i4];
		U1I = x[i12];
		U2R = x[i2];
		U2I = x[i10];
		U3R = x[i6];
		U3I = x[i14];
		U5R = x[i5];
		U5I = x[i13];
		U6R = x[i3];
		U6I = x[i11];
		U7R = x[i7];
		U7I = x[i15];
		T1R = U1R;
		T2R = U2R;
		T3R = U3R;
		T5R = U5R;
		T6R = U6R;
		T7R = U7R;
		U1R = T1R * C1 - U1I * S1;
		U3R = T3R * C3 - U3I * S3;
		U5R = T5R * C5 - U5I * S5;
		U7R = T7R * C7 - U7I * S7;
		U1I = T1R * S1 + U1I * C1;
		U3I = T3R * S3 + U3I * C3;
		U5I = T5R * S5 + U5I * C5;
		U7I = T7R * S7 + U7I * C7;
		U2R = T2R * C2 - U2I * S2;
		U6R = T6R * C6 - U6I * S6;
		U2I = T2R * S2 + U2I * C2;
		U6I = T6R * S6 + U6I * C6;
		T0R = U0R;
		T0I = U0I;
		T2R = U4R;
		T2I = -U4I;
		T1R = U2R + U6R;
		T1I = U2I + U6I;
		T3R = U2R - U6R;
		T3I = U2I - U6I;
		T4R = U1R + U5R;
		T4I = U1I + U5I;
		T6R = U1R - U5R;
		T6I = U1I - U5I;
		T5R = U3R + U7R;
		T5I = U3I + U7I;
		T7R = U3R - U7R;
		T7I = U3I - U7I;
		U0R = T0R + T1R;
		U0I = T0I + T1I;
		U2R = T2R + T3I;
		U2I = T2I - T3R;
		U4R = T0R - T1R;
		U4I = T0I - T1I;
		U6R = T2R - T3I;
		U6I = T2I + T3R;
		U1R = T4R + T5R;
		U1I = T4I + T5I;
		U3R = T6R + T7I;
		U3I = T6I - T7R;
		U5R = T4R - T5R;
		U5I = T4I - T5I;
		U7R = T6R - T7I;
		U7I = T6I + T7R;
		T3R = U3R;
		T5R = U5R;
		T7R = U7R;
		U3R = (T3R + U3I) * M_SQRT1_2;
		U3I = (-T3R + U3I) * M_SQRT1_2;
		U5R = U5I;
		U5I = -T5R;
		U7R = (-T7R + U7I) * M_SQRT1_2;
		U7I = (-T7R - U7I) * M_SQRT1_2;
		T0R = U0R;
		T0I = U0I;
		T2R = U2R;
		T2I = U2I;
		T4R = U4R;
		T4I = U4I;
		T6R = U6R;
		T6I = U6I;
		U0R = T0R + U1R;
		U0I = T0I + U1I;
		U1R = T0R - U1R;
		U1I = T0I - U1I;
		U2R = T2R + U3R;
		U2I = T2I + U3I;
		U3R = T2R - U3R;
		U3I = T2I - U3I;
		U4R = T4R + U5R;
		U4I = T4I + U5I;
		U5R = T4R - U5R;
		U5I = T4I - U5I;
		U6R = T6R + U7R;
		U6I = T6I + U7I;
		U7R = T6R - U7R;
		U7I = T6I - U7I;
		x[i0] = U0R;
		x[i15] = U0I;
		x[i1] = U2R;
		x[i14] = U2I;
		x[i2] = U4R;
		x[i13] = U4I;
		x[i3] = U6R;
		x[i12] = U6I;
		x[i11] = U1R;
		x[i4] = -U1I;
		x[i10] = U3R;
		x[i5] = -U3I;
		x[i9] = U5R;
		x[i6] = -U5I;
		x[i8] = U7R;
		x[i7] = -U7I;
	}
}

template<class T>
void fft_dit_real_split2(T* x, int N) {
	constexpr int N1 = 8;
	const int N2 = N / N1;
	T T0R, T0I, T1R, T1I, T2R, T2I, T3R, T3I;
	T T4R, T4I, T5R, T5I, T6R, T6I, T7R, T7I;
	T U0R, U0I, U1R, U1I, U2R, U2I, U3R, U3I;
	T U4R, U4I, U5R, U5I, U6R, U6I, U7R, U7I;
	T C1, C2, C3, C5, C6, C7, S1, S2, S3, S5, S6, S7;
	int i0, i1, i2, i3, i4, i5, i6, i7;
	int i8, i9, i10, i11, i12, i13, i14, i15;
	int j1, j2, j3, j5, j6, j7;

	if (N <= 1) {
		return;
	} else if (N == 2) {
		T0R = x[0];
		T1R = x[1];
		U0R = T0R + T1R;
		U1R = T0R - T1R;
		x[0] = U0R;
		x[1] = U1R;
		return;
	} else if (N == 4) {
		U0R = x[0];
		U1R = x[2];
		U2R = x[1];
		U3R = x[3];
		T0R = U0R + U2R;
		T2R = U0R - U2R;
		T1R = U1R + U3R;
		T3R = U1R - U3R;
		U0R = T0R + T1R;
		U1R = T2R;
		U1I = -T3R;
		U2R = T0R - T1R;
		x[0] = U0R;
		x[1] = U1R;
		x[2] = U2R;
		x[3] = U1I;
		return;
	}
	fft_dit_real_split2(x, 4 * N2);
	for (int n1 = 4; n1 < N1; n1++) {
		fft_dit_real_split2(x + n1 * N2, N2);
	}
	const auto& W = twiddles(N);
	i0 = 0;
	i1 = i0 + N2;
	i2 = i1 + N2;
	i3 = i2 + N2;
	i4 = i3 + N2;
	i5 = i4 + N2;
	i6 = i5 + N2;
	i7 = i6 + N2;
	U0R = x[i0];
	U1R = x[i4];
	U2R = x[i2];
	U3R = x[i6];
	U4R = x[i1];
	U5R = x[i5];
	U6R = x[i3];
	U7R = x[i7];
	T0R = U0R;
	T2R = U4R;
	T1R = U2R + U6R;
	T3R = U2R - U6R;
	T4R = U1R + U5R;
	T6R = U1R - U5R;
	T5R = U3R + U7R;
	T7R = U3R - U7R;
	U0R = T0R + T1R;
	U2R = T2R;
	U2I = -T3R;
	U4R = T0R - T1R;
	U1R = T4R + T5R;
	U3R = T6R;
	U3I = -T7R;
	U5R = T4R - T5R;
	T3R = U3R;
	T5R = U5R;
	U3R = (T3R + U3I) * M_SQRT1_2;
	U3I = (-T3R + U3I) * M_SQRT1_2;
	U5I = -T5R;
	T0R = U0R;
	T2R = U2R;
	T2I = U2I;
	T4R = U4R;
	T4I = U4I;
	T6R = U2R;
	T6I = -U2I;
	U0R = T0R + U1R;
	U1R = T0R - U1R;
	U2R = T2R + U3R;
	U2I = T2I + U3I;
	U4R = T4R;
	U4I = U5I;
	U6R = T6R - U3R;
	U6I = T6I + U3I;
	x[i0] = U0R;
	x[i1] = U2R;
	x[i2] = U4R;
	x[i3] = U6R;
	x[i4] = U1R;
	x[i5] = U6I;
	x[i6] = U4I;
	x[i7] = U2I;
	if (N2 / 2) {
		constexpr double c1 = cos((1.0 / 8.0) * M_PI);
		constexpr double c3 = cos((3.0 / 8.0) * M_PI);
		i0 = N2 / 2;
		i1 = i0 + N2;
		i2 = i1 + N2;
		i3 = i2 + N2;
		i4 = i3 + N2;
		i5 = i4 + N2;
		i6 = i5 + N2;
		i7 = i6 + N2;
		U0R = x[i0];
		U1R = x[i4];
		U2R = x[i2];
		U3R = x[i6];
		U4R = x[i1];
		U5R = x[i5];
		U6R = x[i3];
		U7R = x[i7];
		T1R = U1R;
		T2R = U2R;
		T3R = U3R;
		T5R = U5R;
		T6R = U6R;
		T7R = U7R;
		U1R = T1R * c1;
		U3R = T3R * c3;
		U5R = T5R * c3;
		U7R = T7R * c1;
		U1I = T1R * c3;
		U3I = T3R * c1;
		U5I = T5R * c1;
		U7I = T7R * c3;
		T0R = U0R;
		T0I = -U4R;
		T1R = (U2R - U6R) * M_SQRT1_2;
		T3R = (U2R + U6R) * M_SQRT1_2;
		T4R = U1R - U5R;
		T4I = -U1I - U5I;
		T6R = U1R + U5R;
		T6I = -U1I + U5I;
		T5R = U3R - U7R;
		T5I = -U3I - U7I;
		T7R = U3R + U7R;
		T7I = -U3I + U7I;
		U0R = T0R + T1R;
		U0I = -T0I - T3R;
		U2R = T0R - T1R;
		U2I = T0I - T3R;
		U6R = T0R + T1R;
		U6I = T0I + T3R;
		U1R = T4R + T5R;
		U1I = T4I + T5I;
		U3R = T6R + T7I;
		U3I = T6I - T7R;
		U5R = T4R - T5R;
		U5I = T4I - T5I;
		U7R = T6R - T7I;
		U7I = T6I + T7R;
		T3R = U3R;
		T5R = U5R;
		T7R = U7R;
		U3R = (T3R + U3I) * M_SQRT1_2;
		U3I = (-T3R + U3I) * M_SQRT1_2;
		U5R = U5I;
		U5I = -T5R;
		U7R = (-T7R + U7I) * M_SQRT1_2;
		U7I = (-T7R - U7I) * M_SQRT1_2;
		U0R = U0R + U1R;
		U0I = U0I + U1I;
		U4R = U2R + U5R;
		U4I = -U2I + U5I;
		U2R = U2R + U3R;
		U2I = U2I + U3I;
		U6R = U6R + U7R;
		U6I = U6I + U7I;
		x[i0] = U0R;
		x[i7] = U0I;
		x[i1] = U2R;
		x[i6] = U2I;
		x[i2] = U4R;
		x[i5] = U4I;
		x[i3] = U6R;
		x[i4] = U6I;
	}
	for (int k2 = 1; k2 < N2 / 2; k2++) {
		j1 = k2;
		j2 = 2 * k2;
		j3 = 3 * k2;
		j5 = 5 * k2;
		j6 = 6 * k2;
		j7 = 7 * k2;
		C1 = W[j1].real();
		C3 = W[j3].real();
		C5 = W[j5].real();
		C7 = W[j7].real();
		S1 = W[j1].imag();
		S3 = W[j3].imag();
		S5 = W[j5].imag();
		S7 = W[j7].imag();
		i0 = k2;
		i1 = i0 + N2;
		i2 = i1 + N2;
		i3 = i2 + N2;
		i4 = i3 + N2;
		i5 = i4 + N2;
		i6 = i5 + N2;
		i7 = i6 + N2;
		i8 = N2 - k2;
		i9 = i8 + N2;
		i10 = i9 + N2;
		i11 = i10 + N2;
		i12 = i11 + N2;
		i13 = i12 + N2;
		i14 = i13 + N2;
		i15 = i14 + N2;
		U0R = x[i0];
		U4R = x[i8];
		U4I = x[i1];
		U0I = x[i9];
		U1R = x[i4];
		U1I = x[i12];
		U2R = x[i2];
		U2I = x[i10];
		U3R = x[i6];
		U3I = x[i14];
		U5R = x[i5];
		U5I = x[i13];
		U6R = x[i3];
		U6I = x[i11];
		U7R = x[i7];
		U7I = x[i15];
		T1R = U1R;
		T2R = U2R;
		T3R = U3R;
		T5R = U5R;
		T6R = U6R;
		T7R = U7R;
		U1R = T1R * C1 - U1I * S1;
		U3R = T3R * C3 - U3I * S3;
		U5R = T5R * C5 - U5I * S5;
		U7R = T7R * C7 - U7I * S7;
		U1I = T1R * S1 + U1I * C1;
		U3I = T3R * S3 + U3I * C3;
		U5I = T5R * S5 + U5I * C5;
		U7I = T7R * S7 + U7I * C7;
		T0R = U0R;
		T0I = U0I;
		T2R = U4R;
		T2I = -U4I;
		T1R = U2R + U6R;
		T1I = U2I + U6I;
		T3R = U2R - U6R;
		T3I = U2I - U6I;
		T4R = U1R + U5R;
		T4I = U1I + U5I;
		T6R = U1R - U5R;
		T6I = U1I - U5I;
		T5R = U3R + U7R;
		T5I = U3I + U7I;
		T7R = U3R - U7R;
		T7I = U3I - U7I;
		U0R = T0R + T1R;
		U0I = T0I + T1I;
		U2R = T2R + T3I;
		U2I = T2I - T3R;
		U4R = T0R - T1R;
		U4I = T0I - T1I;
		U6R = T2R - T3I;
		U6I = T2I + T3R;
		U1R = T4R + T5R;
		U1I = T4I + T5I;
		U3R = T6R + T7I;
		U3I = T6I - T7R;
		U5R = T4R - T5R;
		U5I = T4I - T5I;
		U7R = T6R - T7I;
		U7I = T6I + T7R;
		T3R = U3R;
		T5R = U5R;
		T7R = U7R;
		U3R = (T3R + U3I) * M_SQRT1_2;
		U3I = (-T3R + U3I) * M_SQRT1_2;
		U5R = U5I;
		U5I = -T5R;
		U7R = (-T7R + U7I) * M_SQRT1_2;
		U7I = (-T7R - U7I) * M_SQRT1_2;
		T0R = U0R;
		T0I = U0I;
		T2R = U2R;
		T2I = U2I;
		T4R = U4R;
		T4I = U4I;
		T6R = U6R;
		T6I = U6I;
		U0R = T0R + U1R;
		U0I = T0I + U1I;
		U1R = T0R - U1R;
		U1I = T0I - U1I;
		U2R = T2R + U3R;
		U2I = T2I + U3I;
		U3R = T2R - U3R;
		U3I = T2I - U3I;
		U4R = T4R + U5R;
		U4I = T4I + U5I;
		U5R = T4R - U5R;
		U5I = T4I - U5I;
		U6R = T6R + U7R;
		U6I = T6I + U7I;
		U7R = T6R - U7R;
		U7I = T6I - U7I;
		x[i0] = U0R;
		x[i15] = U0I;
		x[i1] = U2R;
		x[i14] = U2I;
		x[i2] = U4R;
		x[i13] = U4I;
		x[i3] = U6R;
		x[i12] = U6I;
		x[i11] = U1R;
		x[i4] = -U1I;
		x[i10] = U3R;
		x[i5] = -U3I;
		x[i9] = U5R;
		x[i6] = -U5I;
		x[i8] = U7R;
		x[i7] = -U7I;
	}
}

void fft_2pow(double* x, int N) {
	scramble(x, N);
	fft_dit_real_split1((fft_simd4*) x, N / 4);
}

