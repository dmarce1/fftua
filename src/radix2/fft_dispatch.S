#include      "common.h"


              .global        fft_dispatch

              .text


fft_dispatch: bsf            N, %rax
              bsf            N0, %rbx
              sub            %rax, %rbx
              cmp            $0, %rbx
              je             r2
lbl:          test           $1, %rax
              jnz            r8
              jmp            r4
r2:           call           radix2
              ret
r4:           call           radix4
              ret
r8:           call           radix8
              ret
