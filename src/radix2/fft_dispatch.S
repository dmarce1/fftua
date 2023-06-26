#include      "common.h"


              .global        fft_dispatch

              .text


fft_dispatch: bsf            N, %rax
              test           $1, %rax
              jnz            rodd

              call           radix4
              ret

rodd:         call           radix2
              ret
