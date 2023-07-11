#include      "ymm.h"
#define       X              %r12
#define       Y              %r11
#define       NHI            %r10
#define       nhi            %r9
#define       k2             %r8

              .global        dit_rn_iter_n16


              .text

dit_rn_iter_n16:
              push           %r12
              mov            %rdi, X
              mov            %rsi, Y
              mov            %rdx, NHI
              xor            nhi, nhi
              vmovapd        TWO, two
nhi_loop16:   imul           $16, nhi, %rax
              lea            (X, %rax, 8), %rsi
              lea            (Y, %rax, 8), %rdi
              vmovapd        (%rsi), er0
              vmovapd        32(%rsi), er1
              vmovapd        64(%rsi), er2
              vmovapd        96(%rsi), er3
              vmovapd        (%rdi), ei0
              vmovapd        32(%rdi), ei1
              vmovapd        64(%rdi), ei2
              vmovapd        96(%rdi), ei3
              call           transpose_lo
              call           butterfly4
              call           transpose_lo
              vmovapd        TWO, two
              vmovapd        tw16_1r, cos1
              vmovapd        tw16_1i, sin1
              vmovapd        tw16_2r, cos2
              vmovapd        tw16_2i, sin2
              call           tw_butterfly4
              vmovapd        er0, (%rsi)
              vmovapd        er3, 32(%rsi)
              vmovapd        er1, 64(%rsi)
              vmovapd        er2, 96(%rsi)
              vmovapd        ei0, (%rdi)
              vmovapd        ei3, 32(%rdi)
              vmovapd        ei1, 64(%rdi)
              vmovapd        ei2, 96(%rdi)
              inc            nhi
              cmp            nhi, NHI
              jne            nhi_loop16
              pop            %r12
              ret
tw16_1r:      .double        +1.000000000000000000
              .double        +0.923879532511286756
              .double        +0.707106781186547524
              .double        +0.382683432365089771
tw16_1i:      .double        -0.000000000000000000
              .double        -0.382683432365089771
              .double        -0.707106781186547524
              .double        -0.923879532511286756
tw16_2r:      .double        +1.000000000000000000
              .double        +0.707106781186547524
              .double        +0.000000000000000000
              .double        -0.707106781186547524
tw16_2i:      .double        +0.000000000000000000
              .double        -0.707106781186547524
              .double        -1.000000000000000000
              .double        -0.707106781186547524
TWO:          .double        2.0
              .double        2.0
              .double        2.0
              .double        2.0
