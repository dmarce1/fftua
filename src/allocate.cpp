#include <memory>
#include <cassert>

void * operator new(std::size_t n) {
	void* memptr = nullptr;
	if (posix_memalign(&memptr, 32, n) != 0) {
		assert(false);
		abort();
	}
	return memptr;
}

void operator delete(void * p) {
	free(p);
}

void *operator new[](std::size_t s) {
	void* memptr = nullptr;
	if (posix_memalign(&memptr, 32, s) != 0) {
		assert(false);
		abort();
	}
	return memptr;
}

void operator delete[](void *p) {
	free(p);
}
