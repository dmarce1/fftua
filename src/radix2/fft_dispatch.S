#include      "common.h"


              .global        fft_dispatch

              .text


fft_dispatch:
r4:           call           radix2
              ret
