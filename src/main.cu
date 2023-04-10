#include "gpufft.hpp"

#include <stdio.h>
#include <math.h>
#include <unordered_map>

struct fft_plan {
	int N;
	complex<float>* W;
	complex<float>* X;
	complex<float>* Y;
};

__global__ void fft_radix2_ct(fft_plan plan);

__global__ void fft_radix2(fft_plan plan) {
	const int& tid = threadIdx.x;
	const int& bsize = blockDim.x;
	const int& N = plan.N;
	const auto* const w = plan.W;
	auto* x0 = plan.X;
	auto* x1 = plan.Y;
	const int No2 = N >> 1;
	int s = No2;
	int s2 = N;
	while (((long) s2 * (long) s2) > N) {
		for (int i = tid; i < s; i += bsize) {
			int k = i;
			int j = i;
			for (int m = 0; m < No2; m += s) {
				const auto& z0 = x0[k];
				const auto z1 = x0[k + s] * w[m];
				x1[j] = z0 + z1;
				x1[j + No2] = z0 - z1;
				k += s2;
				j += s;
			}
		}
		s >>= 1;
		s2 >>= 1;
		auto* const tmp = x0;
		x0 = x1;
		x1 = tmp;
		__syncthreads();
	}
	while (s >= 1) {
		const int q = s2 * tid;
		const int a = s * tid;
		const int b = s * bsize;
		const int c = s2 * bsize;
		for (int i = 0; i < s; i++) {
			int k = i + q;
			int j = i + a;
			for (int m = s * tid; m < No2; m += b) {
				const auto& z0 = x0[k];
				const auto z1 = x0[k + s] * w[m];
				x1[j] = z0 + z1;
				x1[j + No2] = z0 - z1;
				k += c;
				j += b;
			}
		}
		s >>= 1;
		s2 >>= 1;
		auto* const tmp = x0;
		x0 = x1;
		x1 = tmp;
		__syncthreads();
	}
	if (plan.X != x0) {
		for (int k = tid; k < N; k += bsize) {
			plan.X[k] = x0[k];
		}
	}
	__syncthreads();
}

fft_plan fft_create_plan(int N) {
	fft_plan plan;
	plan.N = N;
	CUDA_CHECK(cudaMallocManaged(&plan.W, sizeof(complex<float> ) * N));
	CUDA_CHECK(cudaMallocManaged(&plan.X, sizeof(complex<float> ) * N));
	CUDA_CHECK(cudaMallocManaged(&plan.Y, sizeof(complex<float> ) * N));
	if (plan.W == nullptr) {
		printf("Memory allocation failed %s %i\n", __FILE__, __LINE__);
	}
	if (plan.X == nullptr) {
		printf("Memory allocation failed %s %i\n", __FILE__, __LINE__);
	}
	if (plan.Y == nullptr) {
		printf("Memory allocation failed %s %i\n", __FILE__, __LINE__);
	}
	for (int n = 0; n < N; n++) {
		const float theta = -2.0 * M_PI * n / N;
		plan.W[n].real() = cos(theta);
		plan.W[n].imag() = sin(theta);
	}
	return plan;
}

void fft_destroy_plan(fft_plan plan) {
	CUDA_CHECK(cudaFree(plan.W));
	CUDA_CHECK(cudaFree(plan.X));
	CUDA_CHECK(cudaFree(plan.Y));
}

int round_up(int n, int m) {
	return (((n - 1) / m) + 1) * m;
}

void gpufft(std::vector<complex<float>>& x) {
	const int N = x.size();
	const int nthreads = std::min(512, round_up(sqrt(N) * 2, 32));
	static std::unordered_map<int, fft_plan> plans;
	if (plans.find(N) == plans.end()) {
		plans[N] = fft_create_plan(N);
	}
	auto* y = plans[N].X;
	for (int n = 0; n < N; n++) {
		y[n] = x[n];
	}
	fft_radix2<<<1, nthreads>>>(plans[N]);
	CUDA_CHECK(cudaDeviceSynchronize());
	for (int n = 0; n < N; n++) {
		x[n] = y[n];
	}
}

float rand1() {
	return 1.0 - 2.0 * (rand() + 0.5) / (RAND_MAX + 1.0);
}

int main(int argc, char **argv) {
	CUDA_CHECK(cudaDeviceSetLimit(cudaLimitStackSize, 16 * 1024));
	timer tm1;
	timer tm2;
	double t1 = 0.0;
	double t2 = 0.0;
	for (int N = 2; N < 128 * 1024 * 1024; N *= 2) {
		double err1 = 0.0;
		double err2 = 0.0;
		for (int i = 0; i < 8; i++) {
			std::vector<complex<float>> X(N);
			for (int n = 0; n < N; n++) {
				X[n].real() = rand1();
				X[n].imag() = rand1();
			}
			auto Y = X;
			auto X0 = X;
			tm2.start();
			fftw(Y);
			tm2.stop();
			tm1.start();
			gpufft(X);
			tm1.stop();
			//	t1 *= 2.0;
			if( i == 0) {
				tm1.reset();
				tm2.reset();
			}
			t1 += tm1.read();
			t2 += tm2.read();
			for (int n = 0; n < N; n++) {
				err1 = err1 + abs(X[n] - Y[n]);
				//	printf("%i : %e %e | %e %e\n", n, X[n].real(), X[n].imag(), Y[n].real(), Y[n].imag());
			}
		}
		printf("%i %e %e %e %e\n", N, err1 / N/8, t1, t2, t2 / t1);
	}
	return 0;
}
