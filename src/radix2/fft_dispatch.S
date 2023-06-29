#include      "common.h"


              .global        fft_dispatch



              .text


fft_dispatch: lea            (radix), %rax
              bsf            N, %rcx
              mov            (%rax, %rcx, 1), %al
              cmp            $8, %al
              je             r8
              cmp            $4, %al
              je             r4
r2:           call           radix2
              ret
r4:           call           radix4
              ret
r8:           call           radix8
              ret


                            .data

              .align         8

radix:        .byte          2
              .byte          4
              .byte          8
              .byte          2
              .byte          4
              .byte          8
              .byte          2
              .byte          4
              .byte          8
              .byte          2
              .byte          4
              .byte          8
              .byte          2
              .byte          4
              .byte          8
              .byte          2
              .byte          4
              .byte          8
              .byte          2
              .byte          4
              .byte          8
              .byte          2
              .byte          4
              .byte          8
              .byte          2
              .byte          4
              .byte          8
              .byte          2
              .byte          4
              .byte          8
              .byte          2
              .byte          4
              .byte          8
              .byte          2
              .byte          4
              .byte          8
              .byte          2
              .byte          4
              .byte          8
              .byte          2
              .byte          4
              .byte          8
              .byte          2
              .byte          4
              .byte          8
              .byte          2
              .byte          4
              .byte          8
              .byte          2
              .byte          4
              .byte          8
              .byte          2
              .byte          4
              .byte          8
              .byte          2
              .byte          4
              .byte          8
              .byte          2
              .byte          4
              .byte          8
              .byte          2
              .byte          4
              .byte          8

