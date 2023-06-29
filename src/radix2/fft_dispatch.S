#include      "common.h"


              .global        fft_dispatch

              .text


fft_dispatch: cmp            $8, N
              je             r2
              cmp            $16, N
              je             r4
              cmp            $32, N
              je             r8
              cmp            $64, N
              je             r4
              cmpl           $128, N
              je             r4
              jmp            r8
r2:           call           radix2
              ret
r4:           call           radix4
              ret
r8:           call           radix8
              ret
