/*
 * types.hpp
 *
 *  Created on: Apr 6, 2023
 *      Author: dmarce1
 */

#ifndef TYPES333_HPP_
#define TYPES333_HPP_

#include <immintrin.h>

typedef double fft_v2df __attribute__ ((vector_size (16)));
typedef double fft_v4df __attribute__ ((vector_size (32)));

class fft_simd2 {
	fft_v2df a;
public:
	fft_simd2& operator=(const fft_simd2& other) = default;
	inline fft_simd2() = default;
	inline fft_simd2& operator=(double b) {
		a = b - (fft_v2df ) { 0.0, 0.0 };
		return *this;
	}
	inline fft_simd2(double b) {
		a = b - (fft_v2df ) { 0.0, 0.0 };
	}
	inline fft_simd2 operator*(const fft_simd2& b) const {
		fft_simd2 c;
		c.a = a * b.a;
		return c;
	}
	inline fft_simd2 operator+(const fft_simd2& b) const {
		fft_simd2 c;
		c.a = a + b.a;
		return c;
	}
	inline fft_simd2 operator-(const fft_simd2& b) const {
		fft_simd2 c;
		c.a = a - b.a;
		return c;
	}
	inline fft_simd2 operator-() const {
		fft_simd2 c;
		c.a = (fft_v2df ) { 0.0, 0.0 } - a;
		return c;
	}
	inline fft_simd2& operator*=(const fft_simd2& b) {
		*this = *this * b;
		return *this;
	}
	inline fft_simd2& operator+=(const fft_simd2& b) {
		*this = *this + b;
		return *this;
	}
	inline fft_simd2& operator-=(const fft_simd2& b) {
		*this = *this - b;
		return *this;
	}
	inline double operator[](int i) const {
		return a[i];
	}
	inline double& operator[](int i) {
		return a[i];
	}
};

#include <limits>

class fft_simd4 {
	fft_v4df a;
public:
	fft_simd4& operator=(const fft_simd4& other) = default;
	inline fft_simd4() {
		a[0] = std::numeric_limits<double>::signaling_NaN();
		a[1] = std::numeric_limits<double>::signaling_NaN();
		a[2] = std::numeric_limits<double>::signaling_NaN();
		a[3] = std::numeric_limits<double>::signaling_NaN();
	}
	inline fft_simd4& operator=(double b) {
		a = b - (fft_v4df ) { 0.0, 0.0, 0.0, 0.0 };
		return *this;
	}
	inline fft_simd4(double b) {
		*this = b;
	}
	inline fft_simd4 operator*(const fft_simd4& b) const {
		fft_simd4 c;
		c.a = a * b.a;
		return c;
	}
	inline fft_simd4 operator+(const fft_simd4& b) const {
		fft_simd4 c;
		c.a = a + b.a;
		return c;
	}
	inline fft_simd4 operator-(const fft_simd4& b) const {
		fft_simd4 c;
		c.a = a - b.a;
		return c;
	}
	inline fft_simd4 operator-() const {
		fft_simd4 c;
		c = fft_simd4(0.0) - *this;
		return c;
	}
	inline fft_simd4& operator*=(const fft_simd4& b) {
		*this = *this * b;
		return *this;
	}
	inline fft_simd4& operator+=(const fft_simd4& b) {
		*this = *this + b;
		return *this;
	}
	inline void gather(double* base, long long* indices) {
		const __m256i i = *((__m256i *) indices);
		(__m256d &) a = _mm256_i64gather_pd(base, i, 1);
	}
	inline fft_simd4& operator-=(const fft_simd4& b) {
		*this = *this - b;
		return *this;
	}
	inline double operator[](int i) const {
		return a[i];
	}
	inline double& operator[](int i) {
		return a[i];
	}
};

template<class T>
class complex {
	T x;
	T y;
public:
	inline T real() const {
		return x;
	}
	inline T& real() {
		return x;
	}
	inline T imag() const {
		return y;
	}
	inline T& imag() {
		return y;
	}
	complex() {
	}
	inline complex(T a) {
		x = a;
		y = 0.0;
	}
	inline complex conj() const {
		complex C;
		C.x = x;
		C.y = -y;
		return C;
	}
	inline friend complex Ix(complex A) {
		complex C;
		C.x = -A.y;
		C.y = A.x;
		return C;
	}
	complex(const complex&) = default;
	template<class U>
	friend inline complex<T> operator*(const complex<T>& A, const complex<U>& B) {
		complex<T> C;
		C.real() = A.real() * B.real() - A.imag() * B.imag();
		C.imag() = A.real() * B.imag() + A.imag() * B.real();
		return C;
	}
	friend inline complex<T> operator+(const complex<T>& A, const complex<T>& B) {
		complex<T> C;
		C.x = A.x + B.x;
		C.y = A.y + B.y;
		return C;
	}
	friend inline complex<T> operator-(const complex<T>& A, const complex<T>& B) {
		complex<T> C;
		C.x = A.x - B.x;
		C.y = A.y - B.y;
		return C;
	}
	friend inline complex<T> operator-(const complex<T>& A) {
		complex<T> C;
		C.x = -A.x;
		C.y = -A.y;
		return C;
	}
	friend inline complex<T> operator*(const complex<T>& A, T& B) {
		complex<T> C;
		C.x = A.x * B;
		C.y = A.y * B;
		return C;
	}
	friend inline complex<T> operator*(T A, const complex<T>& B) {
		return B * A;
	}
	friend inline T abs(const complex<T>& A) {
		return sqrt(A.x * A.x + A.y * A.y);
	}
	template<class U>
	inline complex<T>& operator*=(const complex<U>& other) {
		*this = *this * other;
		return *this;
	}
	inline complex<T>& operator+=(const complex<T>& other) {
		*this = *this + other;
		return *this;
	}
	inline complex<T>& operator-=(const complex<T>& other) {
		*this = *this - other;
		return *this;
	}
	inline complex<T>& operator*=(T other) {
		*this = *this * other;
		return *this;
	}
}
;

#endif /* TYPES_HPP_ */
