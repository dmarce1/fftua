#include "fft.hpp"
#include "timer.hpp"
#include "util.hpp"
#include <cmath>
#include <unordered_map>

fft_method_t fft_select(int N) {
	static thread_local std::unordered_map<int, fft_method_t> values;
	auto i = values.find(N);
	if (i == values.end()) {
		fft_method_t method;
		const auto pfac = prime_factorization(N);
		if (pfac.size() > 1) {
			assert(false);
			abort();
		}
		std::vector<int> allowed;
		for (int r = 2; r <= SFFT_NMAX; r++) {
			if (N >= r && N % r == 0) {
				allowed.push_back(r);
			}
		}
		std::vector<complex<double>> dummy(N);
		double best_time = 1e99;
		int R = -1;
		for (unsigned i = 0; i < allowed.size(); i++) {
			timer tm;
			for (int k = 0; k < 11; k++) {
				fft_cooley_tukey(allowed[i], dummy.data(), N);
				if (k == 0) {
					tm.start();
				}
			}
			tm.stop();
			if (tm.read() < best_time) {
				best_time = tm.read();
				R = allowed[i];
			}
		}
		method.radix = R;
		method.type = FFT_COOLEY_TUKEY;
		values[N] = method;
		i = values.find(N);
	}
	return i->second;
}
