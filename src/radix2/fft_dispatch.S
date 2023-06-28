#include      "common.h"


              .global        fft_dispatch

              .text


fft_dispatch: bsf            N, %rax
              test           $1, %rax
              jz             r4
              jmp            r2
r4:           call           radix4
              ret
r2:           call           radix2
              ret
