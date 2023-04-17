#include "fft.hpp"
#include <memory>
#include <numeric>
#include <unordered_map>

void fft_indices(int* I, int N) {
	if (N <= SFFT_NMAX) {
		return;
	}
	const auto method = fft_select(N);
	switch (method.type) {
	case FFT_COOLEY_TUKEY:
		fft_cooley_tukey_indices(method.radix, I, N);
		break;
	default:
		assert(false);
		abort();
	}

}

const std::vector<int>& fft_indices(int N) {
	static thread_local std::unordered_map<int, std::shared_ptr<std::vector<int>>>values;
	auto i = values.find(N);
	if (i == values.end()) {
		std::vector<int> I(N);
		std::vector<int> J(N);
		std::iota(I.begin(), I.end(), 0);
		auto method = fft_select(N);
		switch(method.type) {
			case FFT_COOLEY_TUKEY:
			fft_cooley_tukey_indices(method.radix, I.data(), N);
			break;
			default:
			assert(false);
			abort();
		}
		for( int n = 0; n < N; n++) {
			J[I[n]] = n;
		}
		values[N] = std::make_shared<std::vector<int>>(std::move(J));
		i = values.find(N);
	}
	return *i->second;
}
