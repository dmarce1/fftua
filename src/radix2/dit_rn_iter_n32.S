#include      "ymm.h"
#define       X              %r8
#define       Y              %r9
#define       NHI            %r10
#define       nhi            %r11
#define       k2             %r12

              .global        dit_rn_iter_n32


              .text

dit_rn_iter_n32:
              push           %r12
              mov            %rdi, X
              mov            %rsi, Y
              mov            %rdx, NHI
              xor            nhi, nhi
              vmovapd        TWO, two
nhi_loop32:   imul           $16, nhi, %rax
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
              vmovapd        128(%rsi), er0
              vmovapd        160(%rsi), er1
              vmovapd        192(%rsi), er2
              vmovapd        224(%rsi), er3
              vmovapd        128(%rdi), ei0
              vmovapd        160(%rdi), ei1
              vmovapd        192(%rdi), ei2
              vmovapd        224(%rdi), ei3
              call           transpose_lo
              call           butterfly4
              call           transpose_lo
              vmovapd        TWO, two
              vmovapd        tw16_1r, cos1
              vmovapd        tw16_1i, sin1
              vmovapd        tw16_2r, cos2
              vmovapd        tw16_2i, sin2
              call           tw_butterfly4
              vmovapd        er0, 128(%rsi)
              vmovapd        er3, 160(%rsi)
              vmovapd        er1, 192(%rsi)
              vmovapd        er2, 224(%rsi)
              vmovapd        ei0, 128(%rdi)
              vmovapd        ei3, 160(%rdi)
              vmovapd        ei1, 192(%rdi)
              vmovapd        ei2, 224(%rdi)
              vmovapd        (%rsi), er0
              vmovapd        128(%rsi), er2
              vmovapd        (%rdi), ei0
              vmovapd        128(%rdi), ei2
              vmovapd        tw32_1r, cos1
              vmovapd        tw32_1i, sin1
              vmovapd        er0, tr
              vmovapd        ei0, ti
              vfmadd231pd    sin1, ei2, tr
              vfnmadd231pd   sin1, er2, ti
              vfnmadd132pd   cos1, tr, er2
              vfnmadd132pd   cos1, ti, ei2
              vfmsub132pd    two, er2, er0
              vfmsub132pd    two, ei2, ei0
              vmovapd        er0, (%rsi)
              vmovapd        er2, 128(%rsi)
              vmovapd        ei0, (%rdi)
              vmovapd        ei2, 128(%rdi)
              vmovapd        32(%rsi), er0
              vmovapd        160(%rsi), er2
              vmovapd        32(%rdi), ei0
              vmovapd        160(%rdi), ei2
              vmovapd        tw32_2r, cos1
              vmovapd        tw32_2i, sin1
              vmovapd        er0, tr
              vmovapd        ei0, ti
              vfmadd231pd    sin1, ei2, tr
              vfnmadd231pd   sin1, er2, ti
              vfnmadd132pd   cos1, tr, er2
              vfnmadd132pd   cos1, ti, ei2
              vfmsub132pd    two, er2, er0
              vfmsub132pd    two, ei2, ei0
              vmovapd        er0, 32(%rsi)
              vmovapd        er2, 160(%rsi)
              vmovapd        ei0, 32(%rdi)
              vmovapd        ei2, 160(%rdi)
              vmovapd        64(%rsi), er0
              vmovapd        192(%rsi), er2
              vmovapd        64(%rdi), ei0
              vmovapd        192(%rdi), ei2
              vmovapd        tw32_3r, cos1
              vmovapd        tw32_3i, sin1
              vmovapd        er0, tr
              vmovapd        ei0, ti
              vfmadd231pd    sin1, ei2, tr
              vfnmadd231pd   sin1, er2, ti
              vfnmadd132pd   cos1, tr, er2
              vfnmadd132pd   cos1, ti, ei2
              vfmsub132pd    two, er2, er0
              vfmsub132pd    two, ei2, ei0
              vmovapd        er0, 64(%rsi)
              vmovapd        er2, 192(%rsi)
              vmovapd        ei0, 64(%rdi)
              vmovapd        ei2, 192(%rdi)
              vmovapd        96(%rsi), er0
              vmovapd        224(%rsi), er2
              vmovapd        96(%rdi), ei0
              vmovapd        224(%rdi), ei2
              vmovapd        tw32_4r, cos1
              vmovapd        tw32_4i, sin1
              vmovapd        er0, tr
              vmovapd        ei0, ti
              vfmadd231pd    sin1, ei2, tr
              vfnmadd231pd   sin1, er2, ti
              vfnmadd132pd   cos1, tr, er2
              vfnmadd132pd   cos1, ti, ei2
              vfmsub132pd    two, er2, er0
              vfmsub132pd    two, ei2, ei0
              vmovapd        er0, 96(%rsi)
              vmovapd        er2, 224(%rsi)
              vmovapd        ei0, 96(%rdi)
              vmovapd        ei2, 224(%rdi)
              inc            nhi
              cmp            nhi, NHI
              jne            nhi_loop32
              pop            %r12
              ret
TWO:          .double        2.0
              .double        2.0
              .double        2.0
              .double        2.0
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
tw32_1r:      .double        1.000000000000000000
              .double        0.980785280403230000
              .double        0.923879532511287000
              .double        0.831469612302545000
tw32_2r:      .double        0.707106781186548000
              .double        0.555570233019602000
              .double        0.382683432365090000
              .double        0.195090322016128000
tw32_3r:      .double        0.000000000000000000
              .double        -0.195090322016128000
              .double        -0.382683432365090000
              .double        -0.555570233019602000
tw32_4r:      .double        -0.707106781186547000
              .double        -0.831469612302545000
              .double        -0.923879532511287000
              .double        -0.980785280403230000
tw32_1i:      .double        0.000000000000000000
              .double        -0.195090322016128000
              .double        -0.382683432365090000
              .double        -0.555570233019602000
tw32_2i:      .double        -0.707106781186547000
              .double        -0.831469612302545000
              .double        -0.923879532511287000
              .double        -0.980785280403230000
tw32_3i:      .double        -1.000000000000000000
              .double        -0.980785280403230000
              .double        -0.923879532511287000
              .double        -0.831469612302545000
tw32_4i:      .double        -0.707106781186548000
              .double        -0.555570233019602000
              .double        -0.382683432365090000
              .double        -0.195090322016129000

