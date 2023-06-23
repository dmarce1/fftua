#include "fftu.hpp"
#include "util.hpp"
#include <memory>

/*
 const auto& W = twiddles((1 << (tbe - tbb)) * N1);
 const int bsz = 1 << bb;
 const int csz = 1 << cb;
 int cmask = 0;
 int tmask = ((1 << tbe) - 1) & ~((1 << tbb) - 1);
 cmask |= (1 << cb);
 for (int n1 = 0; n1 < ilogb(N1); n1++) {
 cmask |= (1 << (bb + n1));
 }
 int n = 0;
 do {
 const int r0 = n;
 const int r1 = r0 + bsz;
 const int r2 = r1 + bsz;
 const int r3 = r2 + bsz;
 const int i0 = n + csz;
 const int i1 = i0 + bsz;
 const int i2 = i1 + bsz;
 const int i3 = i2 + bsz;
 const int ti1 = (n & tmask) >> tbb;
 const int ti2 = 2 * ti1;
 const int ti3 = 3 * ti1;
 T t0, t1, t2, t3, t4, t5, t6;
 t1 = x[r1];
 t2 = x[r2];
 t3 = x[r3];
 x[r1] = t1 * W[ti1].real() - x[i1] * W[ti1].imag();
 x[i1] = t1 * W[ti1].imag() + x[i1] * W[ti1].real();
 x[r2] = t2 * W[ti2].real() - x[i2] * W[ti2].imag();
 x[i2] = t2 * W[ti2].imag() + x[i2] * W[ti2].real();
 x[r3] = t3 * W[ti3].real() - x[i3] * W[ti3].imag();
 x[i3] = t3 * W[ti3].imag() + x[i3] * W[ti3].real();
 t0 = x[r0];
 t1 = x[r1];
 t2 = x[i0];
 t3 = x[i1];
 x[r0] += x[r2];
 x[r1] += x[r3];
 x[i0] += x[i2];
 x[i1] += x[i3];
 x[r2] = t0 - x[r2];
 x[r3] = t1 - x[r3];
 x[i2] = t2 - x[i2];
 x[i3] = t3 - x[i3];
 t0 = x[r0];
 t1 = x[i0];
 t2 = x[r1];
 t3 = x[i1];
 t4 = x[r2];
 t5 = x[i2];
 t6 = x[r3];
 x[r0] += x[r1];
 x[i0] += x[i1];
 x[r1] = t4 + x[i3];
 x[i1] = t5 - t6;
 x[r2] = t0 - t2;
 x[i2] = t1 - t3;
 x[r3] = t4 - x[i3];
 x[i3] = t5 + t6;
 n |= cmask;
 n++;
 n &= ~cmask;
 } while (n < N);
 */

/*

 template<int N1, class T>
 void apply_butterfly(T* x, int N, int bb, int cb, int tbb, int tbe) {
 //	printf( "N = %i bb = %i cb = %i tbb = %i tbe = %i \n", N, bb, cb, tbb, tbe);
 constexpr int NC = 2;
 const auto& W = twiddles((1 << (tbe - tbb)) * N1);
 std::array<T, 2 * N1> u;
 const int bsz = 1 << bb;
 const int csz = 1 << cb;
 int cmask = 0;
 int tmask = ((1 << tbe) - 1) & ~((1 << tbb) - 1);
 cmask |= (1 << cb);
 for (int n1 = 0; n1 < ilogb(N1); n1++) {
 cmask |= (1 << (bb + n1));
 }
 int n = 0;
 do {
 const int twi = (n & tmask) >> tbb;
 for (int bi = 0; bi < N1; bi++) {
 for (int ci = 0; ci < NC; ci++) {
 u[2 * bi + ci] = x[n + bi * bsz + ci * csz];
 }
 }
 for (int bi = 1; bi < N1; bi++) {
 const auto& t = W[bi * twi];
 auto& re = u[2 * bi + 0];
 auto& im = u[2 * bi + 1];
 auto tmp = re;
 re = tmp * t.real() - im * t.imag();
 im = tmp * t.imag() + im * t.real();
 }
 sfft_complex<N1>(u.data());
 for (int bi = 0; bi < N1; bi++) {
 for (int ci = 0; ci < NC; ci++) {
 x[n + bi * bsz + ci * csz] = u[2 * bi + ci];
 }
 }
 n |= cmask;
 n++;
 n &= ~cmask;
 } while (n < N);
 }

 template<class T>
 void apply_butterfly2(T* x, int N, int bb, int cb, int tbb, int tbe) {
 //	printf( "N = %i bb = %i cb = %i tbb = %i tbe = %i \n", N, bb, cb, tbb, tbe);
 constexpr int N1 = 2;
 const auto& W = twiddles((1 << (tbe - tbb)) * N1);
 const int bsz = 1 << bb;
 const int csz = 1 << cb;
 int cmask = 0;
 int tmask = ((1 << tbe) - 1) & ~((1 << tbb) - 1);
 cmask |= (1 << cb);
 for (int n1 = 0; n1 < ilogb(N1); n1++) {
 cmask |= (1 << (bb + n1));
 }
 int n = 0;
 do {
 const int twi = (n & tmask) >> tbb;
 const int r0 = n;
 const int i0 = n + csz;
 const int r1 = r0 + bsz;
 const int i1 = i0 + bsz;
 const auto& w = W[twi];
 const auto c = w.real();
 const auto s = w.imag();
 auto xr0 = x[r0];
 auto xi0 = x[i0];
 auto xr1 = x[r1];
 auto xi1 = x[i1];
 auto t = xr1;
 xr1 = t * c - xi1 * s;
 xi1 = t * s + xi1 * c;
 x[r0] += xr1;
 x[i0] += xi1;
 x[r1] = xr0 - xr1;
 x[i1] = xi0 - xi1;
 n |= cmask;
 n++;
 n &= ~cmask;
 } while (n < N);
 }

 template<class T>
 void apply_butterfly4(T* x, int N, int bb, int cb, int tbb, int tbe) {
 //	printf( "N = %i bb = %i cb = %i tbb = %i tbe = %i \n", N, bb, cb, tbb, tbe);
 constexpr int N1 = 4;
 const auto& W = twiddles((1 << (tbe - tbb)) * N1);
 const int bsz = 1 << bb;
 const int csz = 1 << cb;
 int cmask = 0;
 int tmask = ((1 << tbe) - 1) & ~((1 << tbb) - 1);
 cmask |= (1 << cb);
 for (int n1 = 0; n1 < ilogb(N1); n1++) {
 cmask |= (1 << (bb + n1));
 }
 int n = 0;
 do {
 const int r0 = n;
 const int r1 = r0 + bsz;
 const int r2 = r1 + bsz;
 const int r3 = r2 + bsz;
 const int i0 = n + csz;
 const int i1 = i0 + bsz;
 const int i2 = i1 + bsz;
 const int i3 = i2 + bsz;
 const int ti1 = (n & tmask) >> tbb;
 auto xr0 = x[r0];
 auto xr1 = x[r1];
 auto xr2 = x[r2];
 auto xr3 = x[r3];
 auto xi0 = x[i0];
 auto xi1 = x[i1];
 auto xi2 = x[i2];
 auto xi3 = x[i3];
 T t0, t1, t2, t3;
 T c1, c2, c3;
 T s1, s2, s3;
 c1 = W[ti1].real();
 s1 = W[ti1].imag();
 c2 = c1 * c1 - s1 * s1;
 s2 = 2.0 * c1 * s1;
 c3 = c1 * c2 - s1 * s2;
 s3 = c1 * s2 + s1 * c2;
 t1 = xr1;
 t2 = xr2;
 t3 = xr3;
 xr1 = t1 * c1 - xi1 * s1;
 xi1 = t1 * s1 + xi1 * c1;
 xr2 = t2 * c2 - xi2 * s2;
 xi2 = t2 * s2 + xi2 * c2;
 xr3 = t3 * c3 - xi3 * s3;
 xi3 = t3 * s3 + xi3 * c3;
 t0 = xr0;
 t1 = xr1;
 t2 = xi0;
 t3 = xi1;
 xr0 += xr2;
 xr1 += xr3;
 xi0 += xi2;
 xi1 += xi3;
 xr2 = t0 - xr2;
 xr3 = t1 - xr3;
 xi2 = t2 - xi2;
 xi3 = t3 - xi3;
 x[r0] = xr0 + xr1;
 x[i0] = xi0 + xi1;
 x[r1] = xr2 + xi3;
 x[i1] = xi2 - xr3;
 x[r2] = xr0 - xr1;
 x[i2] = xi0 - xi1;
 x[r3] = xr2 - xi3;
 x[i3] = xi2 + xr3;
 n |= cmask;
 n++;
 n &= ~cmask;
 } while (n < N);
 }


 template<class T>
 void apply_transpose(T* x, int N, int tb1, int tb2, int tnb) {
 int cmask = N - 1;
 int tsz = 1 << tnb;
 int d1 = 1 << tb1;
 int d2 = 1 << tb2;
 cmask = 0;
 for (int n = 0; n < tnb; n++) {
 cmask |= (1 << (tb1 + n));
 cmask |= (1 << (tb2 + n));
 }
 int n = 0;
 do {
 for (int n1 = 0; n1 < tsz; n1++) {
 for (int n2 = n1 + 1; n2 < tsz; n2++) {
 const int i = n + n1 * d1 + n2 * d2;
 const int j = n + n2 * d1 + n1 * d2;
 std::swap(x[i], x[j]);
 }
 }
 n |= cmask;
 n++;
 n &= ~cmask;
 } while (n < N);
 }

 template<class T>
 void apply_butterfly_and_transpose4(T* x, int N, int bb, int cb, int tbb, int tbe, int trb) {
 //	printf( "N = %i bb = %i cb = %i tbb = %i tbe = %i trb = %i\n", N, bb, cb, tbb, tbe, trb);
 constexpr int N1 = 4;
 constexpr int NC = 2;
 const auto& W = twiddles((1 << (tbe - tbb)) * N1);
 std::array<std::array<T, 2 * N1>, N1> u;
 std::array<complex<T>, N1> tw;
 const int bsz = 1 << bb;
 const int tsz = 1 << trb;
 const int csz = 1 << cb;
 int cmask = 0;
 int tmask = ((1 << tbe) - 1) & ~((1 << tbb) - 1);
 cmask |= (1 << cb);
 for (int n1 = 0; n1 < ilogb(N1); n1++) {
 cmask |= (1 << (bb + n1));
 cmask |= (1 << (trb + n1));
 }
 int n = 0;
 do {
 for (int tri = 0; tri < N1; tri++) {
 for (int bi = 0; bi < N1; bi++) {
 for (int ci = 0; ci < NC; ci++) {
 const int i = n + bi * bsz + tri * tsz + ci * csz;
 u[tri][2 * bi + ci] = x[i];
 }
 }
 }
 const int ti = (n & tmask) >> tbb;
 const T c1 = W[ti].real();
 const T s1 = W[ti].imag();
 const T c2 = c1 * c1 - s1 * s1;
 const T s2 = 2.0 * c1 * s1;
 const T c3 = c1 * c2 - s1 * s2;
 const T s3 = c1 * s2 + s1 * c2;
 for (int tri = 0; tri < N1; tri++) {
 auto& u0 = u[tri];
 auto xr0 = u0[0];
 auto xi0 = u0[1];
 auto xr1 = u0[2];
 auto xi1 = u0[3];
 auto xr2 = u0[4];
 auto xi2 = u0[5];
 auto xr3 = u0[6];
 auto xi3 = u0[7];
 T t0, t1, t2, t3;
 t1 = xr1;
 t2 = xr2;
 t3 = xr3;
 xr1 = t1 * c1 - xi1 * s1;
 xi1 = t1 * s1 + xi1 * c1;
 xr2 = t2 * c2 - xi2 * s2;
 xi2 = t2 * s2 + xi2 * c2;
 xr3 = t3 * c3 - xi3 * s3;
 xi3 = t3 * s3 + xi3 * c3;
 t0 = xr0;
 t1 = xr1;
 t2 = xi0;
 t3 = xi1;
 xr0 += xr2;
 xr1 += xr3;
 xi0 += xi2;
 xi1 += xi3;
 xr2 = t0 - xr2;
 xr3 = t1 - xr3;
 xi2 = t2 - xi2;
 xi3 = t3 - xi3;
 const int i0 = n + tri * bsz;
 x[i0] = xr0 + xr1;
 x[i0 + csz] = xi0 + xi1;
 x[i0 + tsz] = xr2 + xi3;
 x[i0 + tsz + csz] = xi2 - xr3;
 x[i0 + 2 * tsz] = xr0 - xr1;
 x[i0 + 2 * tsz + csz] = xi0 - xi1;
 x[i0 + 3 * tsz] = xr2 - xi3;
 x[i0 + 3 * tsz + csz] = xi2 + xr3;
 }
 n |= cmask;
 n++;
 n &= ~cmask;
 } while (n < N);
 }

 template<class T>
 void apply_butterfly_and_transpose4_no_twiddle(T* x, int N, int bb, int cb, int trb) {
 //	printf( "N = %i bb = %i cb = %i tbb = %i tbe = %i trb = %i\n", N, bb, cb, tbb, tbe, trb);
 constexpr int N1 = 4;
 constexpr int NC = 2;
 std::array<std::array<T, 2 * N1>, N1> u;
 const int bsz = 1 << bb;
 const int tsz = 1 << trb;
 const int csz = 1 << cb;
 int cmask = 0;
 cmask |= (1 << cb);
 for (int n1 = 0; n1 < ilogb(N1); n1++) {
 cmask |= (1 << (bb + n1));
 cmask |= (1 << (trb + n1));
 }
 int n = 0;
 do {
 for (int tri = 0; tri < N1; tri++) {
 for (int bi = 0; bi < N1; bi++) {
 for (int ci = 0; ci < NC; ci++) {
 const int i = n + bi * bsz + tri * tsz + ci * csz;
 u[tri][2 * bi + ci] = x[i];
 }
 }
 }
 for (int tri = 0; tri < N1; tri++) {
 auto& u0 = u[tri];
 auto xr0 = u0[0];
 auto xi0 = u0[1];
 auto xr1 = u0[2];
 auto xi1 = u0[3];
 auto xr2 = u0[4];
 auto xi2 = u0[5];
 auto xr3 = u0[6];
 auto xi3 = u0[7];
 T t0, t1, t2, t3;
 t0 = xr0;
 t1 = xr1;
 t2 = xi0;
 t3 = xi1;
 xr0 += xr2;
 xr1 += xr3;
 xi0 += xi2;
 xi1 += xi3;
 xr2 = t0 - xr2;
 xr3 = t1 - xr3;
 xi2 = t2 - xi2;
 xi3 = t3 - xi3;
 const int i0 = n + tri * bsz;
 x[i0] = xr0 + xr1;
 x[i0 + csz] = xi0 + xi1;
 x[i0 + tsz] = xr2 + xi3;
 x[i0 + tsz + csz] = xi2 - xr3;
 x[i0 + 2 * tsz] = xr0 - xr1;
 x[i0 + 2 * tsz + csz] = xi0 - xi1;
 x[i0 + 3 * tsz] = xr2 - xi3;
 x[i0 + 3 * tsz + csz] = xi2 + xr3;
 }
 n |= cmask;
 n++;
 n &= ~cmask;
 } while (n < N);
 }

 */

#include <cstring>

template<class T, int N>
class write_buffer {
	std::array<T, N> buf;
	T* last_store;
	T* ptr;
public:
	write_buffer() {
		ptr = nullptr;
	}
	void flush() {
		if (ptr) {
			int cnt = last_store - ptr + 1;
			if (cnt == N) {
				std::memcpy(ptr, buf.data(), N * sizeof(T));
			} else {
				for (int i = 0; i < cnt; i++) {
					ptr[i] = buf[i];
				}
			}
			ptr = nullptr;
		}
	}
	void store(T* addr, T val) {
		if (ptr == nullptr) {
			ptr = addr;
			last_store = addr;
			buf[0] = val;
		} else {
			if (addr - last_store != 1 || last_store - ptr == N - 1) {
				flush();
				store(addr, val);
			} else {
				buf[addr - ptr] = val;
				last_store = addr;
			}
		}
	}
};

inline int mask_inc(int index, int mask) {
	index |= ~mask;
	index++;
	index &= mask;
	return index;
}

inline int mask_neg(int index, int mask) {
	return mask_inc(~(index & mask), mask) | (index & ~mask);
}

#include <cstring>

template<class T>
void shuffle(T* xin, T* xout, int N, int NLO) {
	const int NHI = N / NLO;
	for (int ihi = 0; ihi < NHI; ihi++) {
		for (int ilo = 0; ilo < NLO; ilo++) {
			const int i = ilo + NLO * ihi;
			const int j = ihi + NHI * ilo;
			xout[j] = xin[i];
		}
	}
}

template<class T>
void apply_transpose(T* x, int N, int lobit, int hibit, int tnb) {
	int n1 = 1 << tnb;
	const int& n2 = n1;
	int d1 = 1 << lobit;
	int d2 = 1 << hibit;
	int nhi = N / d2 / n1;
	int nmid = d2 / d1 / n1;
	int nlo = d1;
	constexpr int nchunk = 32;
	if (nlo < nchunk) {
		for (int ihi = 0; ihi < nhi; ihi++) {
			for (int i1 = 0; i1 < n1; i1++) {
				for (int imid = 0; imid < nmid; imid++) {
					for (int i2 = i1 + 1; i2 < n2; i2++) {
						for (int ilo = 0; ilo < nlo; ilo++) {
							const int i = ilo + nlo * (i1 + n1 * (imid + nmid * (i2 + n2 * ihi)));
							const int j = ilo + nlo * (i2 + n1 * (imid + nmid * (i1 + n2 * ihi)));
							std::swap(x[i], x[j]);
						}
					}
				}
			}
		}
	} else {
		std::array<T, nchunk> buffer;
		for (int ihi = 0; ihi < nhi; ihi++) {
			for (int i1 = 0; i1 < n1; i1++) {
				for (int imid = 0; imid < nmid; imid++) {
					for (int i2 = i1 + 1; i2 < n2; i2++) {
						for (int ilo = 0; ilo < nlo; ilo += nchunk) {
							const int i = ilo + nlo * (i1 + n1 * (imid + nmid * (i2 + n2 * ihi)));
							const int j = ilo + nlo * (i2 + n1 * (imid + nmid * (i1 + n2 * ihi)));
							std::memcpy(buffer.begin(), x + i, sizeof(T) * nchunk);
							std::memcpy(x + i, x + j, sizeof(T) * nchunk);
							std::memcpy(x + j, buffer.begin(), sizeof(T) * nchunk);
						}
					}
				}
			}
		}
	}
}

template<int N1, class T>
void apply_first_butterfly_real(T* x, int N, int bb) {
	//printf("N = %i bb = %i tbb = %i tbe = %i \n", N, bb, tbb, tbe);
	std::array<fft_simd4, 2 * N1> u;
	int dn1 = 1 << (bb - ilogb(N1));
	int Mmask = 0;
	for (int n1 = 0; n1 < ilogb(N1); n1++) {
		Mmask |= (1 << (bb + n1));
	}
	Mmask = ~Mmask;
	Mmask >>= ilogb(N1);
	N >>= ilogb(N1);
	int M = 0;
	auto* z = (fft_simd4*) x;
	do {
		for (int n1 = 0; n1 < N1; n1++) {
			const int i = M + n1 * dn1;
			u[n1] = z[i];
		}
		sfft_real<N1>(u.data());
		int i0 = M;
		z[i0] = u[0];
		for (int n1 = 1; n1 < N1 / 2; n1++) {
			const int i0 = M + n1 * dn1;
			const int i1 = M + (N1 - n1) * dn1;
			z[i0] = u[n1];
			z[i1] = u[N1 - n1];
		}
		i0 = M + (N1 / 2) * dn1;
		z[i0] = u[N1 / 2];
		M = mask_inc(M, Mmask);
	} while (M < N);
}

template<int N1, class T>
void apply_butterfly_real(T* x, int N, int bb, int tbb, int tbe) {
	const int N2 = (1 << (tbe - tbb));
	const auto& W = twiddles(N1 * N2);
	std::array<T, 2 * N1> u;
	const int dn1 = 1 << bb;
	int k2mask = ((1 << tbe) - 1) & ~((1 << tbb) - 1);
	int Mk2mask = 0;
	for (int n1 = 0; n1 < ilogb(N1); n1++) {
		Mk2mask |= (1 << (bb + n1));
	}
	Mk2mask = ~Mk2mask;
	int Mk2 = 0;
	do {
		const int k2 = (Mk2 & k2mask) >> tbb;
		if (k2 == 0 || k2 == N2 / 2) {
			for (int n1 = 0; n1 < N1; n1++) {
				const int i = Mk2 + n1 * dn1;
				u[n1] = x[i];
			}
			if (k2 == 0) {
				sfft_real<N1>(u.data());
				int i0 = Mk2;
				x[i0] = u[0];
				for (int n1 = 1; n1 < N1 / 2; n1++) {
					const int i0 = Mk2 + n1 * dn1;
					const int i1 = Mk2 + (N1 - n1) * dn1;
					x[i0] = u[n1];
					x[i1] = u[N1 - n1];
				}
				i0 = Mk2 + (N1 / 2) * dn1;
				x[i0] = u[N1 / 2];
			} else {
				sfft_skew<N1>(u.data());
				for (int n1 = 0; n1 < N1 / 2; n1++) {
					const int i0 = Mk2 + n1 * dn1;
					const int i1 = Mk2 + (N1 - n1 - 1) * dn1;
					x[i0] = u[n1];
					x[i1] = u[N1 - n1 - 1];
				}
			}
		} else if (k2 < N2 / 2) {
			for (int n1 = 0; n1 < N1; n1++) {
				const int i0 = Mk2 + n1 * dn1;
				const int i1 = mask_neg(Mk2 + n1 * dn1, k2mask);
				u[2 * n1] = x[i0];
				u[2 * n1 + 1] = x[i1];
			}
			for (int n1 = 1; n1 < N1; n1++) {
				const auto& t = W[n1 * k2];
				auto& re = u[2 * n1 + 0];
				auto& im = u[2 * n1 + 1];
				auto tmp = re;
				re = tmp * t.real() - im * t.imag();
				im = tmp * t.imag() + im * t.real();
			}
			sfft_complex<N1>(u.data());
			for (int n1 = 0; n1 < N1 / 2; n1++) {
				const int i0 = Mk2 + n1 * dn1;
				const int i1 = mask_neg(Mk2 + n1 * dn1, ~Mk2mask | k2mask);
				x[i0] = u[2 * n1];
				x[i1] = u[2 * n1 + 1];
			}
			for (int n1 = N1 / 2; n1 < N1; n1++) {
				const int i1 = Mk2 + n1 * dn1;
				const int i0 = mask_neg(Mk2 + n1 * dn1, ~Mk2mask | k2mask);
				x[i0] = u[2 * n1];
				x[i1] = -u[2 * n1 + 1];
			}
		}
		Mk2 = mask_inc(Mk2, Mk2mask);
	} while (Mk2 < N);
}

template<int N1, class T>
void apply_butterfly_and_transpose_real(T* x, int N, int bb, int tbb, int tbe, int trb) {
//	printf( "N = %i bb = %i tbb = %i tbe = %i trb = %i\n", N, bb,  tbb, tbe, trb);
	const int N2 = (1 << (tbe - tbb));
	const auto& W = twiddles(N1 * N2);
	std::array<std::array<T, 2 * N1>, N1> u;
	const int dn1 = 1 << bb;
	const int dnt = 1 << trb;
	int k2mask = ((1 << tbe) - 1) & ~((1 << tbb) - 1);
	int Mk2mask = 0;
	int trmask = 0;
	for (int n1 = 0; n1 < ilogb(N1); n1++) {
		Mk2mask |= (1 << (bb + n1));
		trmask |= (1 << (trb + n1));
		Mk2mask |= (1 << (trb + n1));
	}
	Mk2mask = ~Mk2mask;
	int Mk2 = 0;
	do {
		const int k2 = (Mk2 & k2mask) >> tbb;
		if (k2 == 0 || k2 == N2 / 2) {
			for (int ti = 0; ti < N1; ti++) {
				for (int n1 = 0; n1 < N1; n1++) {
					const int i = Mk2 + n1 * dn1 + ti * dnt;
					u[ti][n1] = x[i];
				}
			}
			if (k2 == 0) {
				for (int ti = 0; ti < N1; ti++) {
					sfft_real<N1>(u[ti].data());
					int i0 = Mk2 + ti * dn1;
					x[i0] = u[ti][0];
					for (int n1 = 1; n1 < N1 / 2; n1++) {
						const int i0 = Mk2 + ti * dn1 + n1 * dnt;
						const int i1 = Mk2 + ti * dn1 + (N1 - n1) * dnt;
						x[i0] = u[ti][n1];
						x[i1] = u[ti][N1 - n1];
					}
					i0 = Mk2 + (N1 / 2) * dnt + ti * dn1;
					x[i0] = u[ti][N1 / 2];
				}
			} else {
				for (int ti = 0; ti < N1; ti++) {
					sfft_skew<N1>(u[ti].data());
					for (int n1 = 0; n1 < N1 / 2; n1++) {
						const int i0 = Mk2 + ti * dn1 + n1 * dnt;
						const int i1 = Mk2 + ti * dn1 + (N1 - n1 - 1) * dnt;
						x[i0] = u[ti][n1];
						x[i1] = u[ti][N1 - n1 - 1];
					}
				}
			}
		} else if (k2 < N2 / 2) {
			for (int ti = 0; ti < N1; ti++) {
				for (int n1 = 0; n1 < N1; n1++) {
					const int i0 = Mk2 + n1 * dn1 + ti * dnt;
					const int i1 = mask_neg(Mk2 + n1 * dn1 + ti * dnt, k2mask);
					u[ti][2 * n1] = x[i0];
					u[ti][2 * n1 + 1] = x[i1];
				}
			}
			for (int ti = 0; ti < N1; ti++) {
				for (int n1 = 1; n1 < N1; n1++) {
					const auto& t = W[n1 * k2];
					auto& re = u[ti][2 * n1 + 0];
					auto& im = u[ti][2 * n1 + 1];
					auto tmp = re;
					re = tmp * t.real() - im * t.imag();
					im = tmp * t.imag() + im * t.real();
				}
				sfft_complex<N1>(u[ti].data());
				for (int n1 = 0; n1 < N1; n1++) {
					const int i0 = Mk2 + n1 * dnt + ti * dn1;
					const int i1 = mask_neg(Mk2 + ti * dn1 + n1 * dnt, trmask | k2mask);
					if (i0 < i1) {
						x[i0] = u[ti][2 * n1];
						x[i1] = u[ti][2 * n1 + 1];
					} else {
						x[i1] = u[ti][2 * n1];
						x[i0] = -u[ti][2 * n1 + 1];
					}
				}
			}
		}
		Mk2 = mask_inc(Mk2, Mk2mask);
	} while (Mk2 < N);
}

const std::vector<__m256d >& sines1(int N) {
	using entry_type = std::shared_ptr<std::vector<__m256d>>;
	static std::unordered_map<int, entry_type> cache;
	auto iter = cache.find(N);
	if (iter != cache.end()) {
		return *(iter->second);
	} else {
		std::vector<__m256d> W(N / SIMD_SIZE);
		for (int n = 0; n < N; n++) {
			W[n / SIMD_SIZE][n % SIMD_SIZE] = sin(-2.0 * M_PI * n / N);
		}
		cache[N] = std::make_shared<std::vector<__m256d>>(std::move(W));
		return *(cache[N]);
	}
}

const std::vector<__m256d >& sines2(int N) {
	using entry_type = std::shared_ptr<std::vector<__m256d>>;
	static std::unordered_map<int, entry_type> cache;
	auto iter = cache.find(N);
	if (iter != cache.end()) {
		return *(iter->second);
	} else {
		std::vector<__m256d> W(N / SIMD_SIZE);
		for (int n = 0; n < N; n++) {
			W[n / SIMD_SIZE][n % SIMD_SIZE] = sin(-4.0 * M_PI * n / N);
		}
		cache[N] = std::make_shared<std::vector<__m256d>>(std::move(W));
		return *(cache[N]);
	}
}

const std::vector<__m256d >& sines3(int N) {
	using entry_type = std::shared_ptr<std::vector<__m256d>>;
	static std::unordered_map<int, entry_type> cache;
	auto iter = cache.find(N);
	if (iter != cache.end()) {
		return *(iter->second);
	} else {
		std::vector<__m256d> W(N / SIMD_SIZE);
		for (int n = 0; n < N; n++) {
			W[n / SIMD_SIZE][n % SIMD_SIZE] = sin(-6.0 * M_PI * n / N);
		}
		cache[N] = std::make_shared<std::vector<__m256d>>(std::move(W));
		return *(cache[N]);
	}
}

const std::vector<__m256d >& cosines1(int N) {
	using entry_type = std::shared_ptr<std::vector<__m256d>>;
	static std::unordered_map<int, entry_type> cache;
	auto iter = cache.find(N);
	if (iter != cache.end()) {
		return *(iter->second);
	} else {
		std::vector<__m256d> W(N / SIMD_SIZE);
		for (int n = 0; n < N; n++) {
			W[n / SIMD_SIZE][n % SIMD_SIZE] = cos(-2.0 * M_PI * n / N);
		}
		cache[N] = std::make_shared<std::vector<__m256d>>(std::move(W));
		return *(cache[N]);
	}
}

const std::vector<__m256d >& cosines2(int N) {
	using entry_type = std::shared_ptr<std::vector<__m256d>>;
	static std::unordered_map<int, entry_type> cache;
	auto iter = cache.find(N);
	if (iter != cache.end()) {
		return *(iter->second);
	} else {
		std::vector<__m256d> W(N / SIMD_SIZE);
		for (int n = 0; n < N; n++) {
			W[n / SIMD_SIZE][n % SIMD_SIZE] = cos(-4.0 * M_PI * n / N);
		}
		cache[N] = std::make_shared<std::vector<__m256d>>(std::move(W));
		return *(cache[N]);
	}
}

const std::vector<__m256d >& cosines3(int N) {
	using entry_type = std::shared_ptr<std::vector<__m256d>>;
	static std::unordered_map<int, entry_type> cache;
	auto iter = cache.find(N);
	if (iter != cache.end()) {
		return *(iter->second);
	} else {
		std::vector<__m256d> W(N / SIMD_SIZE);
		for (int n = 0; n < N; n++) {
			W[n / SIMD_SIZE][n % SIMD_SIZE] = cos(-6.0 * M_PI * n / N);
		}
		cache[N] = std::make_shared<std::vector<__m256d>>(std::move(W));
		return *(cache[N]);
	}
}

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

template<class T>
void scramble(T* x, int N) {
	static std::vector<T> y;
	y.resize(N);
	int N1m = N - 1;
	int j = 0;
	for (int i = 0; i < N1m; i++) {
		y[i] = x[j];
		int k = N >> 1;
		while (k <= j) {
			j -= k;
			k >>= 1;
		}
		j += k;
	}
	std::memcpy(x, y.data(), N * sizeof(T));
}

/*void fft_inplace_real(double* x, int N) {
 constexpr int N1 = 4;
 static timer tm1, tm2;
 int N2 = 1;
 int M = N / (N1 * N2);
 const auto& W = twiddles(N);
 const auto& cos1 = cosines1(N);
 const auto& cos2 = cosines2(N);
 const auto& cos3 = cosines3(N);
 const auto& sin1 = sines1(N);
 const auto& sin2 = sines2(N);
 const auto& sin3 = sines3(N);
 if (ilogb(N) % 2 == 1) {
 const int No2 = N / 2;
 for (int m = 0; m < No2; m += SIMD_SIZE) {
 double* xi0 = x + m;
 double* xi1 = xi0 + No2;
 __m256d U0R = _mm256_load_pd(xi0);
 __m256d U1R = _mm256_load_pd(xi1);
 const __m256d T0 = U0R;
 U0R = _mm256_add_pd(T0, U1R);
 U1R = _mm256_sub_pd(T0, U1R);
 _mm256_store_pd(xi0, U0R);
 _mm256_store_pd(xi1, U1R);
 }
 N2 = 2;
 M = N / (N1 * N2);
 _mm_sfence();
 }
 const __m256d Z0 = _mm256_set1_pd(0.0);
 while (N2 < N / N1) {
 const int MN1 = M * N1;
 int _2k = 0, _1k = 0;
 for (int k2 = 0; k2 < N2; k2++) {
 if (k2 == 0) {
 for (int m = 0; m < M; m += SIMD_SIZE) {
 __m256d T0R, T1R, T2R, T3R, U1I;
 double* xi0 = x + m;
 double* xi1 = xi0 + M;
 double* xi2 = xi1 + M;
 double* xi3 = xi2 + M;
 __m256d U0R = _mm256_load_pd(xi0);
 __m256d U1R = _mm256_load_pd(xi1);
 __m256d U2R = _mm256_load_pd(xi2);
 __m256d U3R = _mm256_load_pd(xi3);
 T0R = _mm256_add_pd(U0R, U2R);
 T2R = _mm256_sub_pd(U0R, U2R);
 T1R = _mm256_add_pd(U1R, U3R);
 T3R = _mm256_sub_pd(U3R, U1R);
 U0R = _mm256_add_pd(T0R, T1R);
 U1R = T2R;
 U1I = T3R;
 U2R = _mm256_sub_pd(T0R, T1R);
 _mm256_store_pd(xi0, U0R);
 _mm256_store_pd(xi2, U1R);
 _mm256_store_pd(xi3, U1I);
 _mm256_store_pd(xi1, U2R);
 }
 } else if (k2 == N2 / 2) {
 const __m256d C0 = _mm256_set1_pd(M_SQRT1_2);
 for (int m = 0; m < M; m += SIMD_SIZE) {
 __m256d T1, T2, U1I, T2R, T0R;
 double* xi0 = x + m + MN1;
 double* xi1 = xi0 + M;
 double* xi2 = xi1 + M;
 double* xi3 = xi2 + M;
 __m256d U0R = _mm256_load_pd(xi0);
 __m256d U1R = _mm256_load_pd(xi1);
 __m256d U2R = _mm256_load_pd(xi2);
 __m256d U3R = _mm256_load_pd(xi3);
 T1 = _mm256_sub_pd(U1R, U3R);
 T2 = _mm256_add_pd(U1R, U3R);
 T1 = _mm256_mul_pd(C0, T1);
 T2 = _mm256_mul_pd(C0, T2);
 T2R = U2R;
 T0R = U0R;
 U0R = _mm256_add_pd(T0R, T1);
 U1R = _mm256_add_pd(T2R, T2);
 U1R = _mm256_sub_pd(Z0, U1R);
 U2R = _mm256_sub_pd(T0R, T1);
 U1I = _mm256_sub_pd(T2R, T2);
 _mm256_store_pd(xi0, U0R);
 _mm256_store_pd(xi3, U1R);
 _mm256_store_pd(xi2, U2R);
 _mm256_store_pd(xi1, U1I);
 }
 } else if (k2 < N2 / 2) {
 const int n2k = ~_1k & (N2 - 1);
 const int j1 = k2 * M;
 const int j2 = 2 * j1;
 const int j3 = 3 * j1;
 const __m256d C1 = _mm256_set1_pd(W[j1].real());
 const __m256d C2 = _mm256_set1_pd(W[j2].real());
 const __m256d C3 = _mm256_set1_pd(W[j3].real());
 const __m256d S1 = _mm256_set1_pd(W[j1].imag());
 const __m256d S2 = _mm256_set1_pd(W[j2].imag());
 const __m256d S3 = _mm256_set1_pd(W[j3].imag());
 for (int m = 0; m < M; m += SIMD_SIZE) {
 __m256d T0R, T0I, T1R, T1I, T2R, T2I, T3R, T3I;
 double* xi0 = x + m + MN1 * _2k;
 double* xi1 = x + m + MN1 * n2k;
 double* xi2 = xi0 + M;
 double* xi3 = xi1 + M;
 double* xi4 = xi2 + M;
 double* xi5 = xi3 + M;
 double* xi6 = xi4 + M;
 double* xi7 = xi5 + M;
 __m256d U0R = _mm256_load_pd(xi0);
 __m256d U0I = _mm256_load_pd(xi1);
 __m256d U1R = _mm256_load_pd(xi2);
 __m256d U1I = _mm256_load_pd(xi3);
 __m256d U2R = _mm256_load_pd(xi4);
 __m256d U2I = _mm256_load_pd(xi5);
 __m256d U3R = _mm256_load_pd(xi6);
 __m256d U3I = _mm256_load_pd(xi7);
 T1R = U1R;
 T2R = U2R;
 T3R = U3R;
 T1I = _mm256_mul_pd(U1I, S1);
 T2I = _mm256_mul_pd(U2I, S2);
 T3I = _mm256_mul_pd(U3I, S3);
 U1R = _mm256_fmsub_pd(T1R, C1, T1I);
 U2R = _mm256_fmsub_pd(T2R, C2, T2I);
 U3R = _mm256_fmsub_pd(T3R, C3, T3I);
 T1I = _mm256_mul_pd(U1I, C1);
 T2I = _mm256_mul_pd(U2I, C2);
 T3I = _mm256_mul_pd(U3I, C3);
 U1I = _mm256_fmadd_pd(T1R, S1, T1I);
 U2I = _mm256_fmadd_pd(T2R, S2, T2I);
 U3I = _mm256_fmadd_pd(T3R, S3, T3I);
 T0R = _mm256_add_pd(U0R, U2R);
 T2R = _mm256_sub_pd(U0R, U2R);
 T0I = _mm256_add_pd(U0I, U2I);
 T2I = _mm256_sub_pd(U0I, U2I);
 T1R = _mm256_add_pd(U1R, U3R);
 T3R = _mm256_sub_pd(U3R, U1R);
 T1I = _mm256_add_pd(U1I, U3I);
 T3I = _mm256_sub_pd(U1I, U3I);
 U0R = _mm256_add_pd(T0R, T1R);
 U3R = _mm256_sub_pd(T2R, T3I);
 U2I = _mm256_sub_pd(T1I, T0I);
 U1I = _mm256_add_pd(T2I, T3R);
 U1R = _mm256_add_pd(T2R, T3I);
 U2R = _mm256_sub_pd(T0R, T1R);
 U3I = _mm256_sub_pd(T3R, T2I);
 U0I = _mm256_add_pd(T0I, T1I);
 _mm256_store_pd(xi0, U0R);
 _mm256_store_pd(xi1, U3R);
 _mm256_store_pd(xi2, U2I);
 _mm256_store_pd(xi3, U1I);
 _mm256_store_pd(xi4, U1R);
 _mm256_store_pd(xi5, U2R);
 _mm256_store_pd(xi6, U3I);
 _mm256_store_pd(xi7, U0I);
 }
 }
 _1k = _2k;
 if (k2 != N2 - 1) {
 int l = N2 >> 1;
 while (l <= _2k) {
 _2k -= l;
 l >>= 1;
 }
 _2k += l;
 } else {
 _2k = N2 - 1;
 }
 }
 N2 *= N1;
 M = N / (N1 * N2);
 }
 scramble(x, N);
 {
 double T0R, T1R, T2R, T3R;
 constexpr int i0 = 0;
 const int i1 = i0 + N2;
 const int i2 = i1 + N2;
 const int i3 = i2 + N2;
 auto U0R = x[i0];
 auto U1R = x[i2];
 auto U2R = x[i1];
 auto U3R = x[i3];
 T0R = U0R + U2R;
 T2R = U0R - U2R;
 T1R = U1R + U3R;
 T3R = U1R - U3R;
 x[i0] = T0R + T1R;
 x[i1] = T2R;
 x[i3] = -T3R;
 x[i2] = T0R - T1R;
 }
 if (N2 >= N1) {
 double T1, T2;
 const int i0 = N2 / 2;
 const int i1 = i0 + N2;
 const int i2 = i1 + N2;
 const int i3 = i2 + N2;
 auto U0R = x[i0];
 auto U1R = x[i2];
 auto U2R = x[i1];
 auto U3R = x[i3];
 T1 = M_SQRT1_2 * (U1R - U3R);
 T2 = M_SQRT1_2 * (U1R + U3R);
 x[i0] = U0R + T1;
 x[i3] = -U2R - T2;
 x[i1] = U0R - T1;
 x[i2] = -T2 + U2R;
 }
 for (int k2 = 1; k2 < std::min(SIMD_SIZE, N2 / 2); k2++) {
 const int j1 = k2 / SIMD_SIZE;
 const auto C1 = cos1[j1][k2 % SIMD_SIZE];
 const auto C2 = cos2[j1][k2 % SIMD_SIZE];
 const auto C3 = cos3[j1][k2 % SIMD_SIZE];
 const auto S1 = sin1[j1][k2 % SIMD_SIZE];
 const auto S2 = sin2[j1][k2 % SIMD_SIZE];
 const auto S3 = sin3[j1][k2 % SIMD_SIZE];
 double T0R, T0I, T1R, T1I, T2R, T2I, T3R, T3I;
 const int i0 = k2;
 const int i1 = N2 - k2;
 const int i2 = i0 + N2;
 const int i3 = i1 + N2;
 const int i4 = i2 + N2;
 const int i5 = i3 + N2;
 const int i6 = i4 + N2;
 const int i7 = i5 + N2;
 auto U0R = x[i0];
 auto U0I = x[i1];
 auto U1R = x[i4];
 auto U1I = x[i5];
 auto U2R = x[i2];
 auto U2I = x[i3];
 auto U3R = x[i6];
 auto U3I = x[i7];
 T1R = U1R;
 T2R = U2R;
 T3R = U3R;
 U1R = T1R * C1 - U1I * S1;
 U2R = T2R * C2 - U2I * S2;
 U3R = T3R * C3 - U3I * S3;
 U1I = T1R * S1 + U1I * C1;
 U2I = T2R * S2 + U2I * C2;
 U3I = T3R * S3 + U3I * C3;
 T0R = U0R + U2R;
 T2R = U0R - U2R;
 T0I = U0I + U2I;
 T2I = U0I - U2I;
 T1R = U1R + U3R;
 T3R = U1R - U3R;
 T1I = U1I + U3I;
 T3I = U1I - U3I;
 x[i0] = T0R + T1R;
 x[i7] = T0I + T1I;
 x[i2] = T2R + T3I;
 x[i5] = T2I - T3R;
 x[i3] = T0R - T1R;
 x[i4] = -T0I + T1I;
 x[i1] = T2R - T3I;
 x[i6] = -T2I - T3R;
 }
 for (int k2 = SIMD_SIZE; k2 < N2 / 2; k2 += SIMD_SIZE) {
 const int j1 = k2 / SIMD_SIZE;
 const auto C1 = cos1[j1];
 const auto C2 = cos2[j1];
 const auto C3 = cos3[j1];
 const auto S1 = sin1[j1];
 const auto S2 = sin2[j1];
 const auto S3 = sin3[j1];
 constexpr int permute = (3 << 0) | (2 << 2) | (1 << 4) | (0 << 6);
 __m256d T0R, T0I, T1R, T1I, T2R, T2I, T3R, T3I;
 double* xi0 = x + k2;
 double* xi1 = x + N2 - k2 - SIMD_SIZE + 1;
 double* xi2 = xi0 + N2;
 double* xi3 = xi1 + N2;
 double* xi4 = xi2 + N2;
 double* xi5 = xi3 + N2;
 double* xi6 = xi4 + N2;
 double* xi7 = xi5 + N2;
 __m256d U0R = _mm256_load_pd(xi0);
 __m256d U0I = _mm256_loadu_pd(xi1);
 __m256d U1R = _mm256_load_pd(xi4);
 __m256d U1I = _mm256_loadu_pd(xi5);
 __m256d U2R = _mm256_load_pd(xi2);
 __m256d U2I = _mm256_loadu_pd(xi3);
 __m256d U3R = _mm256_load_pd(xi6);
 __m256d U3I = _mm256_loadu_pd(xi7);
 U0I = _mm256_permute4x64_pd(U0I, permute);
 U1I = _mm256_permute4x64_pd(U1I, permute);
 U2I = _mm256_permute4x64_pd(U2I, permute);
 U3I = _mm256_permute4x64_pd(U3I, permute);
 T1R = U1R;
 T2R = U2R;
 T3R = U3R;
 T1I = _mm256_mul_pd(U1I, S1);
 T2I = _mm256_mul_pd(U2I, S2);
 T3I = _mm256_mul_pd(U3I, S3);
 U1R = _mm256_fmsub_pd(T1R, C1, T1I);
 U2R = _mm256_fmsub_pd(T2R, C2, T2I);
 U3R = _mm256_fmsub_pd(T3R, C3, T3I);
 T1I = _mm256_mul_pd(U1I, C1);
 T2I = _mm256_mul_pd(U2I, C2);
 T3I = _mm256_mul_pd(U3I, C3);
 U1I = _mm256_fmadd_pd(T1R, S1, T1I);
 U2I = _mm256_fmadd_pd(T2R, S2, T2I);
 U3I = _mm256_fmadd_pd(T3R, S3, T3I);
 T0R = _mm256_add_pd(U0R, U2R);
 T2R = _mm256_sub_pd(U0R, U2R);
 T0I = _mm256_add_pd(U0I, U2I);
 T2I = _mm256_sub_pd(U0I, U2I);
 T1R = _mm256_add_pd(U1R, U3R);
 T3R = _mm256_sub_pd(U1R, U3R);
 T1I = _mm256_add_pd(U1I, U3I);
 T3I = _mm256_sub_pd(U1I, U3I);
 _mm256_store_pd(xi0, _mm256_add_pd(T0R, T1R));
 _mm256_store_pd(xi2, _mm256_add_pd(T2R, T3I));
 _mm256_store_pd(xi4, _mm256_sub_pd(T1I, T0I));
 _mm256_store_pd(xi6, _mm256_sub_pd(Z0, _mm256_add_pd(T2I, T3R)));
 _mm256_storeu_pd(xi7, _mm256_permute4x64_pd(_mm256_add_pd(T0I, T1I), permute));
 _mm256_storeu_pd(xi5, _mm256_permute4x64_pd(_mm256_sub_pd(T2I, T3R), permute));
 _mm256_storeu_pd(xi3, _mm256_permute4x64_pd(_mm256_sub_pd(T0R, T1R), permute));
 _mm256_storeu_pd(xi1, _mm256_permute4x64_pd(_mm256_sub_pd(T2R, T3I), permute));
 }

 }*/

inline __m256& dble2m256(double* x) {
	return *((__m256*) x);
}

extern "C" {
void fft_selfsort(double*, double*, double*, int);
}

void fft_inplace_real(double* x, int N) {
	constexpr int N1 = 4;
	int NHI, NMID, N2, TWHI, KLO, KHI;
	const auto& c = cos_twiddles(N);
	const auto& s = sin_twiddles(N);
	fft_selfsort(x, (double*) c.data(), (double*) s.data(), N);
	return;

	/*	void butterfly4_and_tranpose(double* X, double* W, int NHI, int NMID, int N2);
	 void butterfly4_finish(double* X, double* W, int N);
	 void butterfly4(double* X, double* W, int NHI, int N2);
	 void transpose(double* X, int N2);
	 void trivial_butterfly(double* X, int N);
	 static timer tm1, tm2, tm3, tm4, tm5, tm6, tm7, tm8;
	 N2 = NHI = 1;
	 TWHI = N / (N1 * N2);
	 tm1.start();
	 trivial_butterfly(x, N);
	 tm1.stop();
	 NHI *= N1;
	 N2 *= N1;
	 NMID = N / (NHI * N2 * N1 * N1);
	 while (NMID) {
	 TWHI = N / (N1 * N2);
	 tm2.start();
	 butterfly4_and_transpose(x, (double*) w.data(), NHI, NMID, N2);
	 tm2.stop();
	 NHI *= N1;
	 N2 *= N1;
	 NMID = N / (NHI * N2 * N1 * N1);
	 }
	 NMID = N / (NHI * N2 * N1 * N1);
	 while (N2 * N1 <= N / N1) {
	 TWHI = N / (N1 * N2);
	 NHI = TWHI;
	 tm3.start();
	 butterfly4(x, (double*) w.data(), NHI, N2);
	 tm3.stop();
	 N2 *= N1;
	 }
	 NMID = N / (N1 * N1);
	 tm4.start();
	 transpose(x, N / (N1 * N1));
	 tm4.stop();
	 TWHI = N / (N1 * N2);
	 NHI = TWHI;
	 tm5.start();
	 butterfly4_finish(x, (double*) w.data(), N);
	 tm5.stop();
	 //	double tot = tm1.read() + tm2.read() + tm3.read() + tm4.read() + tm5.read();
	 //	printf("1: %e 2: %e 3: %e 4: %e 5: %e\n", tm1.read() / tot, tm2.read() / tot, tm3.read() / tot, tm4.read() / tot, tm5.read() / tot);

	 const auto& cos1 = cosines1(N);
	 const auto& cos2 = cosines2(N);
	 const auto& cos3 = cosines3(N);
	 const auto& sin1 = sines1(N);
	 const auto& sin2 = sines2(N);
	 const auto& sin3 = sines3(N);
	 const auto k2hi = [](int k2) {
	 return k2 >> 2;
	 };
	 const auto k2lo = [](int k2) {
	 return k2 & 0x3;
	 };
	 const __m256d Z0 = _mm256_set1_pd(0.0);
	 const __m256d TWO = _mm256_set1_pd(2.0);
	 const __m256d SQRT12 = _mm256_set1_pd(M_SQRT1_2);

	 const auto trivial_butterfly4 = [&]() {
	 tm1.start();
	 int M = N / N1;
	 for (int i0 = 0; i0 < M; i0 += SIMD_SIZE) {
	 double* xi0 = x + i0;
	 double* xi1 = xi0 + M;
	 double* xi2 = xi1 + M;
	 double* xi3 = xi2 + M;
	 __m256d T0R, T1R, T2R, T3R;
	 __m256d U0R = _mm256_load_pd(xi0);
	 __m256d U1R = _mm256_load_pd(xi1);
	 __m256d U2R = _mm256_load_pd(xi2);
	 __m256d U3R = _mm256_load_pd(xi3);
	 T0R = _mm256_add_pd(U0R, U2R);
	 T2R = _mm256_sub_pd(U0R, U2R);
	 T1R = _mm256_add_pd(U1R, U3R);
	 T3R = _mm256_sub_pd(U3R, U1R);
	 U0R = _mm256_add_pd(T0R, T1R);
	 U2R = _mm256_sub_pd(T0R, T1R);
	 _mm256_store_pd(xi0, U0R);
	 _mm256_store_pd(xi1, T2R);
	 _mm256_store_pd(xi2, U2R);
	 _mm256_store_pd(xi3, T3R);
	 }
	 tm1.stop();
	 };

	 const auto butterfly_and_transpose4 = [&]() {
	 tm2.start();
	 const int NHIN1 = NHI / N1;
	 const int N2N1 = N2 / N1;
	 const int N1MID = N1 * NMID;
	 const int N1N2MID = N2 * N1MID;
	 std::array<std::array<__m256d, 2 * N1>, N1> u;
	 for (int ihi = 0; ihi < NHIN1; ihi++) {
	 for (int imid = 0; imid < NMID; imid++) {
	 int i0 = N2 * (N1 * (imid + NMID * N1 * ihi));
	 for (int n0 = 0; n0 < N1; n0++) {
	 for (int n1 = 0; n1 < N1; n1++) {
	 const double* xi0 = x + N2 * (n0 + N1 * NMID * n1) + i0;
	 u[n0][n1] = _mm256_load_pd(xi0);
	 }
	 }
	 i0 = N2 * (N1 * (imid + NMID * N1 * ihi));
	 for (int n0 = 0; n0 < N1; n0++) {
	 __m256d T0, T1, T2, T3;
	 double* xi0 = x + N2 * N1 * NMID * n0 + i0;
	 double* xi1 = xi0 + N2;
	 double* xi2 = xi1 + N2;
	 double* xi3 = xi2 + N2;
	 auto& u0 = u[n0];
	 __m256d& U0 = u0[0];
	 __m256d& U1 = u0[1];
	 __m256d& U2 = u0[2];
	 __m256d& U3 = u0[3];
	 T0 = _mm256_add_pd(U0, U2);
	 T2 = _mm256_sub_pd(U0, U2);
	 T1 = _mm256_add_pd(U1, U3);
	 T3 = _mm256_sub_pd(U3, U1);
	 U0 = _mm256_add_pd(T0, T1);
	 U2 = _mm256_sub_pd(T0, T1);
	 _mm256_store_pd(xi0, U0);
	 _mm256_store_pd(xi1, T2);
	 _mm256_store_pd(xi3, T3);
	 _mm256_store_pd(xi2, U2);
	 }
	 }
	 if (N2 >= N1) {
	 const int khi0 = k2hi(N2/2);
	 const int klo0 = k2lo(N2/2);
	 for (int imid = 0; imid < NMID; imid++) {

	 imid + NMIDN1 * ihi

	 int i0 = N1 * (khi0 + N2N1 * (N1 * (imid + NMID * (N1 * (ihi + NHIN1 * klo0)))));;
	 for (int n0 = 0; n0 < N1; n0++) {
	 for (int n1 = 0; n1 < N1; n1++) {
	 const double* xi0 = x + N2 * (n0 + N1 * NMID * n1) + i0;
	 u[n0][n1] = _mm256_load_pd(xi0);
	 }
	 }
	 i0 = N1 * (khi0 + N2N1 * (N1 * (imid + NMID * N1 * (ihi + NHIN1 * klo0))));
	 for (int n0 = 0; n0 < N1; n0++) {
	 __m256d T1, T2, T0, T3;
	 double* xi0 = x + N2 * N1 * NMID * n0 + i0;
	 double* xi1 = xi0 + N2;
	 double* xi2 = xi1 + N2;
	 double* xi3 = xi2 + N2;
	 auto& u0 = u[n0];
	 auto& U0 = u0[0];
	 auto& U1 = u0[1];
	 auto& U2 = u0[2];
	 auto& U3 = u0[3];
	 T0 = _mm256_add_pd(U1, U3);
	 T2 = _mm256_sub_pd(U1, U3);
	 T1 = _mm256_mul_pd(SQRT12, T2);
	 T3 = _mm256_mul_pd(SQRT12, T0);
	 T0 = U0;
	 T2 = U2;
	 U0 = _mm256_add_pd(T0, T1);
	 U3 = _mm256_add_pd(T2, T3);
	 U1 = _mm256_sub_pd(T0, T1);
	 U2 = _mm256_sub_pd(T2, T3);
	 U3 = _mm256_sub_pd(Z0, U3);
	 _mm256_store_pd(xi0, U0);
	 _mm256_store_pd(xi1, U1);
	 _mm256_store_pd(xi2, U2);
	 _mm256_store_pd(xi3, U3);
	 }
	 }
	 }
	 for (int imid = 0; imid < NMID; imid++) {
	 for (int k2 = 1; k2 < N2 / 2; k2++) {
	 const int j1 = k2 * TWHI;
	 const int j2 = 2 * j1;
	 const int j3 = 3 * j1;
	 const __m256d C1 = _mm256_set1_pd(w[j1].real());
	 const __m256d C2 = _mm256_set1_pd(w[j2].real());
	 const __m256d C3 = _mm256_set1_pd(w[j3].real());
	 const __m256d S1 = _mm256_set1_pd(w[j1].imag());
	 const __m256d S2 = _mm256_set1_pd(w[j2].imag());
	 const __m256d S3 = _mm256_set1_pd(w[j3].imag());
	 const int khi0 = k2hi(k2);
	 const int klo0 = k2lo(k2);
	 const int khi1 = k2hi(N2 - k2);
	 const int klo1 = k2lo(N2 - k2);
	 const int i0 = N1 * (khi0 + N2N1 * (N1 * (imid + NMID * (N1 * (ihi + NHIN1 * klo0)))));
	 const int i1 = N1 * (khi1 + N2N1 * (N1 * (imid + NMID * (N1 * (ihi + NHIN1 * klo1)))));
	 for (int n0 = 0; n0 < N1; n0++) {
	 for (int n1 = 0; n1 < N1; n1++) {
	 const int di = N2 * (n0 + N1MID * n1);
	 const double* xi0 = x + i0 + di;
	 const double* xi1 = x + i1 + di;
	 u[n0][2 * n1 + 0] = _mm256_load_pd(xi0);
	 u[n0][2 * n1 + 1] = _mm256_load_pd(xi1);
	 }
	 }
	 for (int n0 = 0; n0 < N1; n0++) {
	 __m256d T0R, T0I, T1R, T1I, T2R, T2I, T3R, T3I;
	 const int di = N1N2MID * n0;
	 double* xi0 = x + i0 + di;
	 double* xi4 = x + i1 + di;
	 double* xi1 = xi0 + N2;
	 double* xi2 = xi1 + N2;
	 double* xi3 = xi2 + N2;
	 double* xi5 = xi4 + N2;
	 double* xi6 = xi5 + N2;
	 double* xi7 = xi6 + N2;
	 auto& u0 = u[n0];
	 __m256d& U0R = u0[0];
	 __m256d& U0I = u0[1];
	 __m256d& U1R = u0[2];
	 __m256d& U1I = u0[3];
	 __m256d& U2R = u0[4];
	 __m256d& U2I = u0[5];
	 __m256d& U3R = u0[6];
	 __m256d& U3I = u0[7];
	 T1R = U1R;
	 T0R = _mm256_fmsub_pd(S2, U2I, U0R);
	 T0R = _mm256_fmsub_pd(C2, U2R, T0R);
	 T0I = _mm256_fmadd_pd(C2, U2I, U0I);
	 T0I = _mm256_fmadd_pd(S2, U2R, T0I);
	 T2R = _mm256_fmsub_pd(TWO, U0R, T0R);
	 T2I = _mm256_fmsub_pd(TWO, U0I, T0I);
	 T1I = _mm256_mul_pd(U1I, S1);
	 U1R = _mm256_fmsub_pd(T1R, C1, T1I);
	 T1I = _mm256_mul_pd(U1I, C1);
	 U1I = _mm256_fmadd_pd(T1R, S1, T1I);
	 T1R = _mm256_fmsub_pd(S3, U3I, U1R);
	 T1R = _mm256_fmsub_pd(C3, U3R, T1R);
	 T1I = _mm256_fmadd_pd(C3, U3I, U1I);
	 T1I = _mm256_fmadd_pd(S3, U3R, T1I);
	 T3R = _mm256_fmsub_pd(TWO, U1R, T1R);
	 T3I = _mm256_fmsub_pd(TWO, U1I, T1I);
	 T3R = _mm256_sub_pd(Z0, T3R);
	 U0R = _mm256_add_pd(T0R, T1R);
	 U0I = _mm256_add_pd(T0I, T1I);
	 U1R = _mm256_add_pd(T2R, T3I);
	 U1I = _mm256_add_pd(T2I, T3R);
	 U2R = _mm256_sub_pd(T0R, T1R);
	 U2I = _mm256_sub_pd(T1I, T0I);
	 U3R = _mm256_sub_pd(T2R, T3I);
	 U3I = _mm256_sub_pd(T3R, T2I);
	 _mm256_store_pd(xi0, U0R);
	 _mm256_store_pd(xi1, U1R);
	 _mm256_store_pd(xi2, U2I);
	 _mm256_store_pd(xi3, U3I);
	 _mm256_store_pd(xi4, U3R);
	 _mm256_store_pd(xi5, U2R);
	 _mm256_store_pd(xi6, U1I);
	 _mm256_store_pd(xi7, U0I);
	 }
	 }
	 }
	 }
	 tm2.stop();
	 };

	 const auto butterfly_and_transpose2 = [&]() {
	 tm3.start();
	 const int NHIN1 = NHI / N1;
	 const int N2N1 = N2 / N1;
	 std::array<std::array<double, N1>, N1 / 2> u;
	 for (int ihi = 0; ihi < NHIN1; ihi++) {
	 for (int imid = 0; imid < NMID; imid++) {
	 for(int ilo = 0; ilo < N1; ilo++) {
	 const int khi0 = k2hi(0);
	 const int klo0 = k2lo(0);
	 for (int n0 = 0; n0 < N1 / 2; n0++) {
	 for (int n1 = 0; n1 < N1 / 2; n1++) {
	 const int i0 = ilo + N1 * (khi0 + N2N1 * (n0 + (N1/2) * (imid + NMID * (n1 + (N1/2) * (ihi + NHIN1 * klo0)))));
	 u[n0][n1] = x[i0];
	 }
	 }
	 for (int n0 = 0; n0 < (N1 / 2); n0++) {
	 const int i0 = ilo + N1 * (khi0 + N2N1 * ((N1/2) * (imid + NMID * (n0 + (N1/2) * (ihi + NHIN1 * klo0)))));
	 const int i1 = i0 + N2;
	 const auto& u0 = u[n0];
	 auto U0R = u0[0];
	 auto U1R = u0[1];
	 x[i0] = U0R + U1R;
	 x[i1] = U0R - U1R;
	 }
	 }
	 }
	 if (N2 >= (N1 / 2)) {
	 const int khi0 = k2hi(N2/2);
	 const int klo0 = k2lo(N2/2);
	 for (int imid = 0; imid < NMID; imid++) {
	 for(int ilo = 0; ilo < N1; ilo++) {
	 for (int n0 = 0; n0 < (N1 / 2); n0++) {
	 for (int n1 = 0; n1 < (N1 / 2); n1++) {
	 const int i0 = ilo + N1 * (khi0 + N2N1 * (n0 + (N1/2) * (imid + NMID * (n1 + (N1/2) * (ihi + NHIN1 * klo0)))));
	 u[n0][n1] = x[i0];
	 }
	 }
	 for (int n0 = 0; n0 < (N1 / 2); n0++) {
	 const int i0 = ilo + N1 * (khi0 + N2N1 * ((N1/2) * (imid + NMID * (n0 + (N1/2) * (ihi + NHIN1 * klo0)))));
	 const int i1 = i0 + N2;
	 const auto& u0 = u[n0];
	 auto U0R = u0[0];
	 auto U1R = u0[1];
	 x[i0] = U0R;
	 x[i1] = -U1R;
	 }
	 }
	 }
	 }
	 for (int k2 = 1; k2 < N2 / 2; k2++) {
	 const int khi0 = k2hi(k2);
	 const int klo0 = k2lo(k2);
	 const int khi1 = k2hi(N2 - k2);
	 const int klo1 = k2lo(N2 - k2);
	 const int j = k2 * TWHI;
	 const auto C = w[j].real();
	 const auto S = w[j].imag();
	 for (int imid = 0; imid < NMID; imid++) {
	 for(int ilo = 0; ilo < N1; ilo++) {
	 for (int n0 = 0; n0 < (N1 / 2); n0++) {
	 for (int n1 = 0; n1 < (N1 / 2); n1++) {
	 const int i0 = ilo + N1 * (khi0 + N2N1 * (n0 + (N1/2) * (imid + NMID * (n1 + (N1/2) * (ihi + NHIN1 * klo0)))));
	 const int i1 = ilo + N1 * (khi1 + N2N1 * (n0 + (N1/2) * (imid + NMID * (n1 + (N1/2) * (ihi + NHIN1 * klo1)))));
	 u[n0][2 * n1] = x[i0];
	 u[n0][2 * n1 + 1] = x[i1];
	 }
	 }
	 for (int n0 = 0; n0 < (N1 / 2); n0++) {
	 const int i0 = ilo + N1 * (khi0 + N2N1 * ((N1/2) * (imid + NMID * (n0 + (N1/2) * (ihi + NHIN1 * klo0)))));
	 const int i2 = ilo + N1 * (khi1 + N2N1 * ((N1/2) * (imid + NMID * (n0 + (N1/2) * (ihi + NHIN1 * klo1)))));
	 const int i1 = i0 + N2;
	 const int i3 = i2 + N2;
	 const auto& u0 = u[n0];
	 auto U0R = u0[0];
	 auto U0I = u0[1];
	 auto U1R = u0[2];
	 auto U1I = u0[3];
	 auto T = U1R;
	 U1R = T * C - U1I * S;
	 U1I = T * S + U1I * C;
	 x[i0] = U0R + U1R;
	 x[i3] = U0I + U1I;
	 x[i2] = U0R - U1R;
	 x[i1] = -U0I + U1I;
	 }
	 }
	 }
	 }
	 }
	 tm3.stop();
	 };

	 const auto butterfly2 = [&]() {
	 tm4.start();
	 const int NHIN1 = NHI / N1;
	 const int N2N1 = N2 / N1;
	 for (int ihi = 0; ihi < NHIN1; ihi++) {
	 for( int ilo = 0; ilo < N1; ilo++) {
	 const int khi0 = k2hi(0);
	 const int klo0 = k2lo(0);
	 const int i0 = ilo + N1 * (khi0 + N2N1 * ((N1/2) * (ihi + NHIN1 * klo0)));
	 const int i1 = i0 + N2;
	 auto U0R = x[i0];
	 auto U1R = x[i1];
	 x[i0] = U0R + U1R;
	 x[i1] = U0R - U1R;
	 }
	 }
	 for (int ihi = 0; ihi < NHIN1; ihi++) {
	 for( int ilo = 0; ilo < N1; ilo++) {
	 const int khi0 = k2hi(N2/2);
	 const int klo0 = k2lo(N2/2);
	 const int i1 = ilo + N1 * (khi0 + N2N1 * ((N1/2) * (ihi + NHIN1 * klo0))) + N2;
	 x[i1] = -x[i1];
	 }
	 }
	 for (int ihi = 0; ihi < NHIN1; ihi++) {
	 for (int k2 = 1; k2 < N2 / 2; k2++) {
	 const int khi0 = k2hi(k2);
	 const int klo0 = k2lo(k2);
	 const int khi1 = k2hi(N2 - k2);
	 const int klo1 = k2lo(N2 - k2);
	 const int j = k2 * TWHI;
	 const auto C = w[j].real();
	 const auto S = w[j].imag();
	 for( int ilo = 0; ilo < N1; ilo++) {
	 const int i0 = ilo + N1 * (khi0 + N2N1 * ((N1/2) * (ihi + NHIN1 * klo0)));
	 const int i2 = ilo + N1 * (khi1 + N2N1 * ((N1/2) * (ihi + NHIN1 * klo1)));
	 const int i1 = i0 + N2;
	 const int i3 = i2 + N2;
	 auto U0R = x[i0];
	 auto U0I = x[i2];
	 auto U1R = x[i1];
	 auto U1I = x[i3];
	 auto T = U1R;
	 U1R = T * C - U1I * S;
	 U1I = T * S + U1I * C;
	 x[i0] = U0R + U1R;
	 x[i3] = U0I + U1I;
	 x[i2] = U0R - U1R;
	 x[i1] = -U0I + U1I;
	 }
	 }
	 }
	 tm4.stop();
	 };

	 const auto butterfly4 = [&]() {
	 tm5.start();
	 const int NHIN1 = NHI / N1;
	 const int N2N1 = N2 / N1;
	 for (int ihi = 0; ihi < NHIN1; ihi++) {
	 __m256d T0, T1, T2, T3;
	 double* xi0 = x + N2 * N1 * ihi;
	 double* xi1 = xi0 + N2;
	 double* xi2 = xi1 + N2;
	 double* xi3 = xi2 + N2;
	 __m256d U0 = _mm256_load_pd(xi0);
	 __m256d U1 = _mm256_load_pd(xi1);
	 __m256d U2 = _mm256_load_pd(xi2);
	 __m256d U3 = _mm256_load_pd(xi3);
	 T0 = _mm256_add_pd(U0, U2);
	 T2 = _mm256_sub_pd(U0, U2);
	 T1 = _mm256_add_pd(U1, U3);
	 T3 = _mm256_sub_pd(U3, U1);
	 U0 = _mm256_add_pd(T0, T1);
	 U2 = _mm256_sub_pd(T0, T1);
	 _mm256_store_pd(xi0, U0);
	 _mm256_store_pd(xi1, T2);
	 _mm256_store_pd(xi3, T3);
	 _mm256_store_pd(xi2, U2);
	 }
	 if (N2 >= N1) {
	 const int khi0 = k2hi(N2/2);
	 const int klo0 = k2lo(N2/2);
	 for (int ihi = 0; ihi < NHIN1; ihi++) {
	 __m256d T1, T2, T0, T3;
	 double* xi0 = x + N1 * (khi0 + N2N1 * (0 + N1 * (ihi + NHIN1 * klo0)));
	 double* xi1 = xi0 + N2;
	 double* xi2 = xi1 + N2;
	 double* xi3 = xi2 + N2;
	 __m256d U0 = _mm256_load_pd(xi0);
	 __m256d U1 = _mm256_load_pd(xi1);
	 __m256d U2 = _mm256_load_pd(xi2);
	 __m256d U3 = _mm256_load_pd(xi3);
	 T0 = _mm256_add_pd(U1, U3);
	 T2 = _mm256_sub_pd(U1, U3);
	 T1 = _mm256_mul_pd(SQRT12, T2);
	 T3 = _mm256_mul_pd(SQRT12, T0);
	 T0 = U0;
	 T2 = U2;
	 U0 = _mm256_add_pd(T0, T1);
	 U3 = _mm256_add_pd(T2, T3);
	 U1 = _mm256_sub_pd(T0, T1);
	 U2 = _mm256_sub_pd(T2, T3);
	 U3 = _mm256_sub_pd(Z0, U3);
	 _mm256_store_pd(xi0, U0);
	 _mm256_store_pd(xi1, U1);
	 _mm256_store_pd(xi2, U2);
	 _mm256_store_pd(xi3, U3);
	 }
	 }
	 for (int ihi = 0; ihi < NHIN1; ihi++) {
	 for( int klo = 0; klo < N1; klo++) {
	 for (int khi = 0; khi < N2N1; khi++) {
	 const int k2 = (khi << 2) | klo;
	 if( k2 == 0 || k2 >= N2/2) {
	 continue;
	 }
	 const int& khi0 = khi;
	 const int& klo0 = klo;
	 const int khi1 = k2hi(N2 - k2);
	 const int klo1 = k2lo(N2 - k2);
	 const int j1 = k2 * TWHI;
	 const int j2 = 2 * j1;
	 const int j3 = 3 * j1;
	 const __m256d C1 = _mm256_set1_pd(w[j1].real());
	 const __m256d C2 = _mm256_set1_pd(w[j2].real());
	 const __m256d C3 = _mm256_set1_pd(w[j3].real());
	 const __m256d S1 = _mm256_set1_pd(w[j1].imag());
	 const __m256d S2 = _mm256_set1_pd(w[j2].imag());
	 const __m256d S3 = _mm256_set1_pd(w[j3].imag());
	 __m256d T0R, T0I, T1R, T1I, T2R, T2I, T3R, T3I;
	 double* xi0 = x + N1 * (khi0 + N2 * (ihi + NHIN1 * klo0));
	 double* xi4 = x + N1 * (khi1 + N2 * (ihi + NHIN1 * klo1));
	 double* xi1 = xi0 + N2;
	 double* xi2 = xi1 + N2;
	 double* xi3 = xi2 + N2;
	 double* xi5 = xi4 + N2;
	 double* xi6 = xi5 + N2;
	 double* xi7 = xi6 + N2;
	 __m256d U0R = _mm256_load_pd(xi0);
	 __m256d U0I = _mm256_load_pd(xi4);
	 __m256d U1R = _mm256_load_pd(xi1);
	 __m256d U1I = _mm256_load_pd(xi5);
	 __m256d U2R = _mm256_load_pd(xi2);
	 __m256d U2I = _mm256_load_pd(xi6);
	 __m256d U3R = _mm256_load_pd(xi3);
	 __m256d U3I = _mm256_load_pd(xi7);
	 T1R = U1R;
	 T0R = _mm256_fmsub_pd(S2, U2I, U0R);
	 T0R = _mm256_fmsub_pd(C2, U2R, T0R);
	 T0I = _mm256_fmadd_pd(C2, U2I, U0I);
	 T0I = _mm256_fmadd_pd(S2, U2R, T0I);
	 T2R = _mm256_fmsub_pd(TWO, U0R, T0R);
	 T2I = _mm256_fmsub_pd(TWO, U0I, T0I);
	 T1I = _mm256_mul_pd(U1I, S1);
	 U1R = _mm256_fmsub_pd(T1R, C1, T1I);
	 T1I = _mm256_mul_pd(U1I, C1);
	 U1I = _mm256_fmadd_pd(T1R, S1, T1I);
	 T1R = _mm256_fmsub_pd(S3, U3I, U1R);
	 T1R = _mm256_fmsub_pd(C3, U3R, T1R);
	 T1I = _mm256_fmadd_pd(C3, U3I, U1I);
	 T1I = _mm256_fmadd_pd(S3, U3R, T1I);
	 T3R = _mm256_fmsub_pd(TWO, U1R, T1R);
	 T3I = _mm256_fmsub_pd(TWO, U1I, T1I);
	 T3R = _mm256_sub_pd(Z0, T3R);
	 U0R = _mm256_add_pd(T0R, T1R);
	 U0I = _mm256_add_pd(T0I, T1I);
	 U1R = _mm256_add_pd(T2R, T3I);
	 U1I = _mm256_add_pd(T2I, T3R);
	 U2R = _mm256_sub_pd(T0R, T1R);
	 U2I = _mm256_sub_pd(T1I, T0I);
	 U3R = _mm256_sub_pd(T2R, T3I);
	 U3I = _mm256_sub_pd(T3R, T2I);
	 _mm256_store_pd(xi0, U0R);
	 _mm256_store_pd(xi1, U1R);
	 _mm256_store_pd(xi2, U2I);
	 _mm256_store_pd(xi3, U3I);
	 _mm256_store_pd(xi4, U3R);
	 _mm256_store_pd(xi5, U2R);
	 _mm256_store_pd(xi6, U1I);
	 _mm256_store_pd(xi7, U0I);
	 }
	 }
	 }
	 tm5.stop();
	 };

	 const auto butterfly4_final = [&]() {
	 tm6.start();
	 for (int ihi = 0; ihi < NHI; ihi++) {
	 double T0R, T1R, T2R, T3R;
	 const int i0 = N2 * N1 * ihi;
	 const int i1 = i0 + N2;
	 const int i2 = i1 + N2;
	 const int i3 = i2 + N2;
	 auto U0R = x[i0];
	 auto U1R = x[i1];
	 auto U2R = x[i2];
	 auto U3R = x[i3];
	 T0R = U0R + U2R;
	 T2R = U0R - U2R;
	 T1R = U1R + U3R;
	 T3R = U1R - U3R;
	 x[i0] = T0R + T1R;
	 x[i1] = T2R;
	 x[i3] = -T3R;
	 x[i2] = T0R - T1R;
	 }
	 if (N2 >= N1) {
	 for (int ihi = 0; ihi < NHI; ihi++) {
	 double T1, T2;
	 const int i0 = N2 / 2 + N2 * N1 * ihi;
	 const int i1 = i0 + N2;
	 const int i2 = i1 + N2;
	 const int i3 = i2 + N2;
	 auto U0R = x[i0];
	 auto U1R = x[i1];
	 auto U2R = x[i2];
	 auto U3R = x[i3];
	 T1 = M_SQRT1_2 * (U1R - U3R);
	 T2 = M_SQRT1_2 * (U1R + U3R);
	 x[i0] = U0R + T1;
	 x[i3] = -U2R - T2;
	 x[i1] = U0R - T1;
	 x[i2] = -T2 + U2R;
	 }
	 }
	 for (int ihi = 0; ihi < NHI; ihi++) {
	 for (int k2 = 1; k2 < std::min(SIMD_SIZE,N2 / 2); k2++) {
	 double T0R, T0I, T1R, T1I, T2R, T2I, T3R, T3I;
	 const int i0 = k2 + N2 * N1 * ihi;
	 const int i1 = i0 + N2;
	 const int i2 = i1 + N2;
	 const int i3 = i2 + N2;
	 const int i4 = i0 + N2 - 2 * k2;
	 const int i5 = i4 + N2;
	 const int i6 = i5 + N2;
	 const int i7 = i6 + N2;
	 const int j1 = k2 * TWHI;
	 const int j2 = 2 * j1;
	 const int j3 = 3 * j1;
	 auto U0R = x[i0];
	 auto U1R = x[i1];
	 auto U2R = x[i2];
	 auto U3R = x[i3];
	 auto U0I = x[i4];
	 auto U1I = x[i5];
	 auto U2I = x[i6];
	 auto U3I = x[i7];
	 const auto C1 = w[j1].real();
	 const auto C2 = w[j2].real();
	 const auto C3 = w[j3].real();
	 const auto S1 = w[j1].imag();
	 const auto S2 = w[j2].imag();
	 const auto S3 = w[j3].imag();
	 T1R = U1R;
	 T2R = U2R;
	 T3R = U3R;
	 U1R = T1R * C1 - U1I * S1;
	 U2R = T2R * C2 - U2I * S2;
	 U3R = T3R * C3 - U3I * S3;
	 U1I = T1R * S1 + U1I * C1;
	 U2I = T2R * S2 + U2I * C2;
	 U3I = T3R * S3 + U3I * C3;
	 T0R = U0R + U2R;
	 T2R = U0R - U2R;
	 T0I = U0I + U2I;
	 T2I = U0I - U2I;
	 T1R = U1R + U3R;
	 T3R = U1R - U3R;
	 T1I = U1I + U3I;
	 T3I = U1I - U3I;
	 x[i0] = T0R + T1R;
	 x[i7] = T0I + T1I;
	 x[i1] = T2R + T3I;
	 x[i6] = T2I - T3R;
	 x[i5] = T0R - T1R;
	 x[i2] = -T0I + T1I;
	 x[i4] = T2R - T3I;
	 x[i3] = -T2I - T3R;
	 }
	 constexpr unsigned int P = (3 << 0) | (2 << 2) | (1 << 4) | (0 << 6);
	 for (int k2 = SIMD_SIZE; k2 < N2 / 2; k2 += SIMD_SIZE) {
	 __m256d T0R, T0I, T1R, T1I, T2R, T2I, T3R, T3I;
	 double* xi0 = x + k2 + N2 * N1 * ihi;
	 double* xi1 = xi0 + N2;
	 double* xi2 = xi1 + N2;
	 double* xi3 = xi2 + N2;
	 double* xi4 = xi0 + N2 - 2 * k2 - SIMD_SIZE + 1;
	 double* xi5 = xi4 + N2;
	 double* xi6 = xi5 + N2;
	 double* xi7 = xi6 + N2;
	 const int j = k2 * TWHI / SIMD_SIZE;
	 const auto C1 = cos1[j];
	 const auto S1 = sin1[j];
	 const auto C2 = cos2[j];
	 const auto S2 = sin2[j];
	 const auto C3 = cos3[j];
	 const auto S3 = sin3[j];
	 __m256d U0R = _mm256_load_pd(xi0);
	 __m256d U0I = _mm256_loadu_pd(xi4);
	 __m256d U1R = _mm256_load_pd(xi1);
	 __m256d U1I = _mm256_loadu_pd(xi5);
	 __m256d U2R = _mm256_load_pd(xi2);
	 __m256d U2I = _mm256_loadu_pd(xi6);
	 __m256d U3R = _mm256_load_pd(xi3);
	 __m256d U3I = _mm256_loadu_pd(xi7);
	 U0I = _mm256_permute4x64_pd(U0I, P);
	 U1I = _mm256_permute4x64_pd(U1I, P);
	 U2I = _mm256_permute4x64_pd(U2I, P);
	 U3I = _mm256_permute4x64_pd(U3I, P);
	 T1R = U1R;
	 T0R = _mm256_fmsub_pd(S2, U2I, U0R);
	 T0R = _mm256_fmsub_pd(C2, U2R, T0R);
	 T0I = _mm256_fmadd_pd(C2, U2I, U0I);
	 T0I = _mm256_fmadd_pd(S2, U2R, T0I);
	 T2R = _mm256_fmsub_pd(TWO, U0R, T0R);
	 T2I = _mm256_fmsub_pd(TWO, U0I, T0I);
	 T1I = _mm256_mul_pd(U1I, S1);
	 U1R = _mm256_fmsub_pd(T1R, C1, T1I);
	 T1I = _mm256_mul_pd(U1I, C1);
	 U1I = _mm256_fmadd_pd(T1R, S1, T1I);
	 T1R = _mm256_fmsub_pd(S3, U3I, U1R);
	 T1R = _mm256_fmsub_pd(C3, U3R, T1R);
	 T1I = _mm256_fmadd_pd(C3, U3I, U1I);
	 T1I = _mm256_fmadd_pd(S3, U3R, T1I);
	 T3R = _mm256_fmsub_pd(TWO, U1R, T1R);
	 T3I = _mm256_fmsub_pd(TWO, U1I, T1I);
	 T3R = _mm256_sub_pd(Z0, T3R);
	 U0R = _mm256_add_pd(T0R, T1R);
	 U0I = _mm256_add_pd(T0I, T1I);
	 U1R = _mm256_add_pd(T2R, T3I);
	 U1I = _mm256_add_pd(T2I, T3R);
	 U2R = _mm256_sub_pd(T0R, T1R);
	 U2I = _mm256_sub_pd(T1I, T0I);
	 U3R = _mm256_sub_pd(T2R, T3I);
	 U3I = _mm256_sub_pd(T3R, T2I);
	 U0I = _mm256_permute4x64_pd(U0I, P);
	 U1I = _mm256_permute4x64_pd(U1I, P);
	 U2R = _mm256_permute4x64_pd(U2R, P);
	 U3R = _mm256_permute4x64_pd(U3R, P);
	 _mm256_store_pd(xi0, U0R);
	 _mm256_storeu_pd(xi7, U0I);
	 _mm256_store_pd(xi1, U1R);
	 _mm256_storeu_pd(xi6, U1I);
	 _mm256_storeu_pd(xi5, U2R);
	 _mm256_store_pd(xi2, U2I);
	 _mm256_storeu_pd(xi4, U3R);
	 _mm256_store_pd(xi3, U3I);
	 }
	 }
	 tm6.stop();
	 };

	 N2 = NHI = 1;
	 TWHI = N / (N1 * N2);
	 trivial_butterfly4();
	 NHI *= N1;
	 N2 *= N1;
	 NMID = N / (NHI * N2 * N1 * N1);
	 while (NMID) {
	 TWHI = N / (N1 * N2);
	 butterfly_and_transpose4();
	 NHI *= N1;
	 N2 *= N1;
	 NMID = N / (NHI * N2 * N1 * N1);
	 }
	 NMID = N / (NHI * N2 * N1 * N1);
	 if (NHI * N2 == N / 8) {
	 NMID = N1 / 2;
	 TWHI = N / (N1 / 2 * N2);
	 butterfly_and_transpose2();
	 N2 *= (N1 / 2);
	 NHI = N / N2 / (N1 / 2) / (N1 / 2);
	 tm7.start();
	 for (int ihi = 0; ihi < NHI; ihi++) {
	 const int i = N2 * (1 + (N1 / 2) * ((N1 / 2) * ihi));
	 const int j = N2 * (N1 / 2) * (1 + (N1 / 2) * ihi);
	 for (int ilo = 0; ilo < N2; ilo++) {
	 std::swap(x[i + ilo], x[j + ilo]);
	 }
	 }
	 tm7.stop();
	 } else {
	 if (NHI * N2 == N / 2) {
	 TWHI = N / (N1 / 2 * N2);
	 NHI = TWHI;
	 butterfly2();
	 N2 *= N1 / 2;
	 }
	 }
	 while (N2 * N1 <= N / N1) {
	 TWHI = N / (N1 * N2);
	 NHI = TWHI;
	 butterfly4();
	 N2 *= N1;
	 }
	 NMID = N / (N1 * N1);
	 tm8.start();
	 for (int n1 = 0; n1 < N1; n1++) {
	 for (int n2 = n1 + 1; n2 < N1; n2++) {
	 for (int n0 = 0; n0 < NMID; n0++) {
	 const int i = n2 + N1 * (n0 + NMID * n1);
	 const int j = n1 + N1 * (n0 + NMID * n2);
	 std::swap(x[i], x[j]);
	 }
	 }
	 }
	 tm8.stop();
	 TWHI = N / (N1 * N2);
	 NHI = TWHI;
	 butterfly4_final();
	 */
}

template<int N1, class T>
void apply_butterfly(T* x, int N, int bb, int cb, int tbb, int tbe) {
//	printf( "N = %i bb = %i cb = %i tbb = %i tbe = %i \n", N, bb, cb, tbb, tbe);
	constexpr int NC = 2;
	const auto& W = twiddles(N1 * (1 << (tbe - tbb)));
	std::array<T, 2 * N1> u;
	const int bsz = 1 << bb;
	const int csz = 1 << cb;
	int tmask = ((1 << tbe) - 1) & ~((1 << tbb) - 1);
	int cmask = 0;
	cmask |= (1 << cb);
	for (int n1 = 0; n1 < ilogb(N1); n1++) {
		cmask |= (1 << (bb + n1));
	}
	int n = 0;
	do {
		const int twi = (n & tmask) >> tbb;
		for (int bi = 0; bi < N1; bi++) {
			for (int ci = 0; ci < NC; ci++) {
				u[2 * bi + ci] = x[n + bi * bsz + ci * csz];
			}
		}
		for (int bi = 1; bi < N1; bi++) {
			const auto& t = W[bi * twi];
			auto& re = u[2 * bi + 0];
			auto& im = u[2 * bi + 1];
			auto tmp = re;
			re = tmp * t.real() - im * t.imag();
			im = tmp * t.imag() + im * t.real();
		}
		sfft_complex<N1>(u.data());
		for (int bi = 0; bi < N1; bi++) {
			for (int ci = 0; ci < NC; ci++) {
				x[n + bi * bsz + ci * csz] = u[2 * bi + ci];
			}
		}
		n |= cmask;
		n++;
		n &= ~cmask;
	} while (n < N);
}

template<int N1, class T>
void apply_butterfly_and_transpose(T* x, int N, int bb, int cb, int tbb, int tbe, int trb) {
//	printf( "N = %i bb = %i cb = %i tbb = %i tbe = %i trb = %i\n", N, bb, cb, tbb, tbe, trb);
	constexpr int NC = 2;
	const auto& W = twiddles(N1 * (1 << (tbe - tbb)));
	std::array<std::array<T, 2 * N1>, N1> u;
	std::array<complex<T>, N1> tw;
	const int bsz = 1 << bb;
	const int tsz = 1 << trb;
	const int csz = 1 << cb;
	int cmask = 0;
	int tmask = ((1 << tbe) - 1) & ~((1 << tbb) - 1);
	cmask |= (1 << cb);
	for (int n1 = 0; n1 < ilogb(N1); n1++) {
		cmask |= (1 << (bb + n1));
		cmask |= (1 << (trb + n1));
	}
	int n = 0;
	do {
		const int twi = (n & tmask) >> tbb;
		for (int tri = 0; tri < N1; tri++) {
			for (int bi = 0; bi < N1; bi++) {
				for (int ci = 0; ci < NC; ci++) {
					const int i = n + bi * bsz + tri * tsz + ci * csz;
					u[tri][2 * bi + ci] = x[i];
				}
			}
		}
		if (tmask) {
			for (int k1 = 1; k1 < N1; k1++) {
				tw[k1] = W[k1 * twi];
			}
		}
		for (int tri = 0; tri < N1; tri++) {
			if (tmask) {
				for (int bi = 1; bi < N1; bi++) {
					auto& re = u[tri][2 * bi + 0];
					auto& im = u[tri][2 * bi + 1];
					auto tmp = re;
					re = tmp * tw[bi].real() - im * tw[bi].imag();
					im = tmp * tw[bi].imag() + im * tw[bi].real();
				}
			}
			sfft_complex<N1>(u[tri].data());
			for (int bi = 0; bi < N1; bi++) {
				for (int ci = 0; ci < NC; ci++) {
					x[n + tri * bsz + bi * tsz + ci * csz] = u[tri][2 * bi + ci];
				}
			}
		}
		n |= cmask;
		n++;
		n &= ~cmask;
	} while (n < N);

}

void fft_inplace(double* x, int N) {
	if (N <= SFFT_NMAX) {
		sfft_complex(x, N);
		return;
	}
	constexpr int N1 = 4;
	const int highest_bit = ilogb(N) - ilogb(N1) + 1;
	int lobit = 1;
	int hibit = highest_bit;
	while (hibit > lobit + ilogb(N1) - 1) {
		apply_butterfly_and_transpose<N1, double>(x, 2 * N, hibit, 0, 1, lobit, lobit);
		lobit += ilogb(N1);
		hibit -= ilogb(N1);
	}
	if (hibit - lobit == 1) {
		apply_butterfly<2 * N1, double>(x, 2 * N, lobit, 0, 1, lobit);
		lobit += 3;
	}
	if (hibit - lobit == -1) {
		apply_butterfly<N1 / 2, double>(x, 2 * N, lobit, 0, 1, lobit);
		lobit += 1;
	}
	while (lobit <= highest_bit) {
		apply_butterfly<N1, double>(x, 2 * N, lobit, 0, 1, lobit);
		lobit += ilogb(N1);
	}
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
template<int N1, class T>
void fft_inplace_part2(T* X, T* Y, const complex<double>* W, int N, int N2) {
	const int NHI = N / (N1 * N2);
	std::array<T, 2 * N1> u;
	std::array<complex<double>, N1> w;
	for (int nhi = 0; nhi < NHI; nhi++) {
		for (int n2 = 0; n2 < N2; n2++) {
			for (int na = 0; na < N1; na++) {
				const int i = n2 + N2 * (na + N1 * nhi);
				u[2 * na] = X[i];
				u[2 * na + 1] = Y[i];
			}
			for (int na = 1; na < N1; na++) {
				const auto ti = na * n2 * NHI;
				w[na] = W[ti];
			}
			for (int k = 1; k < N1; k++) {
				auto& r = u[2 * k];
				auto& i = u[2 * k + 1];
				auto tmp = r;
				r = tmp * w[k].real() - i * w[k].imag();
				i = tmp * w[k].imag() + i * w[k].real();
			}
			sfft_complex<N1>(u.data());
			for (int na = 0; na < N1; na++) {
				const int i = n2 + N2 * (na + N1 * nhi);
				X[i] = u[2 * na];
				Y[i] = u[2 * na + 1];
			}
		}
	}
}

void fft_inplace1(double* X, double* Y, int N) {
	constexpr int N0 = 2;
	constexpr int N1 = 4;
	constexpr int N1o2 = 2;
	const auto& W = twiddles(N);
	int NHI = 1;
	int N2 = 1;
	int NMID = N / (N1 * N1);
	std::array<std::array<double, 2 * N1>, N1> u;
	std::array<complex<double>, N1> w;

	while (NMID) {
		for (int nhi = 0; nhi < NHI; nhi++) {
			for (int nmid = 0; nmid < NMID; nmid++) {
				for (int n2 = 0; n2 < N2; n2++) {
					for (int na = 0; na < N1; na++) {
						for (int nb = 0; nb < N1; nb++) {
							const int i = n2 + N2 * (na + N1 * (nmid + NMID * (nb + N1 * nhi)));
							u[na][2 * nb] = X[i];
							u[na][2 * nb + 1] = Y[i];
						}
					}
					for (int na = 1; na < N1; na++) {
						const auto ti = na * n2 * NHI * N1 * NMID;
						w[na] = W[ti];
					}
					for (int na = 0; na < N1; na++) {
						for (int k = 1; k < N1; k++) {
							auto& r = u[na][2 * k];
							auto& i = u[na][2 * k + 1];
							auto tmp = r;
							r = tmp * w[k].real() - i * w[k].imag();
							i = tmp * w[k].imag() + i * w[k].real();
						}
						sfft_complex<N1>(u[na].data());
					}
					for (int na = 0; na < N1; na++) {
						for (int nb = 0; nb < N1; nb++) {
							const int i = n2 + N2 * (nb + N1 * (nmid + NMID * (na + N1 * nhi)));
							X[i] = u[na][2 * nb];
							Y[i] = u[na][2 * nb + 1];
						}
					}
				}
			}
		}
		NHI *= N1;
		N2 *= N1;
		NMID /= N1 * N1;
	}
	if (ilogb(N) % 4 == 1) {
		N2 = 1 << (ilogb(N) / 2);
		fft_inplace_part2<N1 / N0, double>(X, Y, W.data(), N, N2);
		N2 *= N1o2;
	} else if (ilogb(N) % 4 == 3) {
		N2 = 1 << (ilogb(N) / 2 - 1);
		fft_inplace_part2<N1 * N0, double>(X, Y, W.data(), N, N2);
		N2 = 1 << (ilogb(N) / 2 + 2);
	}
	while (N2 < N) {
		fft_inplace_part2<N1, double>(X, Y, W.data(), N, N2);
		N2 *= N1;
	}
}

void fft_inplace3(double* X, int N) {
	static std::vector<double> xre, xim;
	xre.resize(N);
	xim.resize(N);
	for (int n = 0; n < N; n++) {
		xre[n] = X[2 * n];
		xim[n] = X[2 * n + 1];
	}
	fft_inplace1(xre.data(), xim.data(), N);
	for (int n = 0; n < N; n++) {
		X[2 * n] = xre[n];
		X[2 * n + 1] = xim[n];
	}
}

void fft_width(double* X, int N) {
	constexpr int N1 = 4;

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

void fft_width2(double* X, int N) {
	constexpr int N1 = 4;
	constexpr int cache_size = 1024;

	int nepoch = ceil(log(N) / log(N1) / floor(log(cache_size) / log(N1)));
	const int nlev = ilogb(N) / ilogb(N1);
	while (nlev % nepoch != 0) {
		nepoch++;
	}
	const int npass = nlev / nepoch;
	int nbutter = cache_size;
	int ngroup = N / nbutter;
	if (N < cache_size) {
		nbutter = N;
		ngroup = 1;
	}
	while (ilogb(ngroup) % ilogb(nbutter) != 0) {
		nbutter /= N1;
		ngroup *= N1;
	}

//	printf("N1 = %i cache_size = %i nlev = %i nepoch = %i npass = %i ngroup = %i nbutter = %i\n", N1, cache_size, nlev, nepoch, npass, ngroup, nbutter);

	scramble_complex(X, N);
	const auto& W = twiddles(N);
	int nghi = ngroup;
	int nglo = 1;
	int digit = 0;
	std::array<double, 2 * N1> u;
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
							for (int bi = 0; bi < N1; bi++) {
								u[2 * bi + 0] = X[2 * (k0 + bi * nblo) + 0];
								u[2 * bi + 1] = X[2 * (k0 + bi * nblo) + 1];
							}
							const int ti = (glo + nglo * blo) * tishft;
							for (int bi = 1; bi < N1; bi++) {
								const auto& t = W[bi * ti];
								auto& re = u[2 * bi + 0];
								auto& im = u[2 * bi + 1];
								auto tmp = re;
								re = tmp * t.real() - im * t.imag();
								im = tmp * t.imag() + im * t.real();
							}
							sfft_complex<N1>(u.data());
							for (int bi = 0; bi < N1; bi++) {
								X[2 * (k0 + bi * nblo) + 0] = u[2 * bi + 0];
								X[2 * (k0 + bi * nblo) + 1] = u[2 * bi + 1];
							}
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
							if (i < j) {
								std::swap(X[2 * i], X[2 * j]);
								std::swap(X[2 * i + 1], X[2 * j + 1]);
							}
						}
					}
				}
			}
		}
		nghi /= nbutter;
		nglo *= nbutter;
		digit += npass;
	}
	int nhi = N / nbutter;
	int nlo = nbutter;
	int jhi = 0;
	for (int ihi = 0; ihi < nhi - 1; ihi++) {
		if (ihi > jhi) {
			for (int n = 0; n < 2 * nbutter; n++) {
				std::swap(X[2 * ihi * nlo + n], X[2 * jhi * nlo + n]);
			}
		}
		int k = nhi >> 1;
		while (k <= jhi) {
			jhi -= k;
			k >>= 1;
		}
		jhi += k;
	}
	for (int ihi = 0; ihi < nhi; ihi++) {
		scramble_complex(X + 2 * ihi * nlo, nbutter);
	}
	scramble_complex(X, N);

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

