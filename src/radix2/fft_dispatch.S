#include      "common.h"


              .global        fft_dispatch

              .text


fft_dispatch:
r4:           call           radix4
              ret
