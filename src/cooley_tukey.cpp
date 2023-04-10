#include "fftu.hpp"
#include "util.hpp"
#include <unordered_map>
#include <cmath>
#include <vector>
#include <stack>
#include <cstring>

void fft_cooley_tukey_indices(int N1, int* I, int N) {
	const int N2 = N / N1;
	std::vector<int> J(N);
	for (int n2 = 0; n2 < N2; n2++) {
		J[n2] = I[N1 * n2];
		J[N / 2 + n2] = I[N1 * n2 + N1 / 2];
		for (int n1 = 1; n1 < N1 / 2; n1++) {
			J[N2 * n1 + n2] = I[N1 * n2 + n1];
		}
		for (int n1 = 1; n1 < N1 / 2; n1++) {
			J[N2 * (N1 - n1) + n2] = I[mod(N1 * n2 - n1, N)];
		}
	}
	for (int n1 = 0; n1 < N1; n1++) {
		fft_indices(J.data() + n1 * N2, N2);
	}
	std::memcpy(I, J.data(), sizeof(int) * N);
}

struct scratch_space {
	std::vector<complex<fft_simd4>> z;
};

static std::stack<scratch_space> spaces;

static void destroy_scratch(scratch_space&& space) {
	spaces.push(std::move(space));
}

static scratch_space create_scratch(int N) {
	scratch_space space;
	if (spaces.size()) {
		space = std::move(spaces.top());
		spaces.pop();
	}
	space.z.resize(N);
	return space;
}

void fft_cooley_tukey(int N1, complex<fft_simd4>* X, int N) {
	const int N2 = N / N1;
	const int No2 = N / 2;
	const auto& W = twiddles(N);
	auto scratch = create_scratch(N1);
	auto& z = scratch.z;
	const int N1o2 = N1 / 2;
	for (int n1 = 0; n1 < N1; n1++) {
		fft(X + N2 * n1, N2);
	}
	for (int k2 = 0; k2 < N2; k2++) {
		z[0] = X[k2];
		z[N1o2] = X[No2 + k2] * W[k2 * N1o2];
		for (int n1 = 1; n1 < N1o2; n1++) {
			const auto& w = W[k2 * n1];
			z[n1] = X[N2 * n1 + k2] * w;
			z[N1 - n1] = X[N2 * (N1 - n1) + k2] * w.conj();
		}
		fft_scramble(z.data(), N1);
		fft(z.data(), N1);
		for (int k1 = 0; k1 < N1; k1++) {
			X[N2 * k1 + k2] = z[k1];
		}
	}
	destroy_scratch(std::move(scratch));
}
