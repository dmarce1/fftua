#include      "ymm.h"
#define       X              %r15
#define       Y              %r14
#define       N2             %r13
#define       Wr             %r12
#define       Wi             %r11
#define       k2             %r10
#define       Wptr           %r9
#define       stack_size     %r8
#define       L1_size        $0x8000


              .global        fft_radix2_complex


              .text

fft_radix2_complex:
              push           %rbx
              push           %r15
              push           %r14
              push           %r13
              push           %r12
              mov            %rdi, X
              mov            %rsi, Y
              mov            %rdx, N2
              shr            $2, N2
              bsf            N2, %rax
              add            $3, %rax
              shl            $4, %rax
              push           %rbp
              mov            %rsp, %rbp
              sub            %rax, %rsp
              mov            %rsp, Wptr
              mov            %rax, stack_size
              bsf            N2, %rdx
              add            $2, %rdx
              mov            $2, %rdi
twiddle_loop: mov            %rdi, %rcx
              push           %rdi
              mov            $1, %rdi
              shl            %rcx, %rdi
              push           %rdx
              push           %r8
              push           %r9
              push           %r10
              push           %r11
              call           twiddles_fwd_complex
              pop            %r11
              pop            %r10
              pop            %r9
              pop            %r8
              pop            %rcx
              pop            %rdi
              mov            %rdi, %rsi
              shl            %rsi
              mov            %rax, (Wptr, %rsi, 8)
              mov            %rdx, 8(Wptr, %rsi, 8)
              mov            %rcx, %rdx
              cmp            %rdi, %rdx
              je             twiddle_esc
              inc            %rdi
              jmp            twiddle_loop
twiddle_esc:  call           next_level
              add            stack_size, %rsp
              mov            %rsp, %rbp
              pop            %rbp
              pop            %r12
              pop            %r13
              pop            %r14
              pop            %r15
              pop            %rbx
              ret
next_level:

cmp            $16, N2
              jl             Nlt64
              push           Y
              push           X
              imul           $3, N2, %rdx
              lea            (X, %rdx, 8), %rax
              lea            (Y, %rdx, 8), %rbx
              push           %rbx
              push           %rax
              imul           $2, N2, %rdx
              lea            (X, %rdx, 8), %rax
              lea            (Y, %rdx, 8), %rbx
              push           %rbx
              push           %rax
              lea            (X, N2, 8), %rax
              lea            (Y, N2, 8), %rbx
              push           %rbx
              push           %rax
              shr            $2, N2
              call           next_level
              pop            X
              pop            Y
              call           next_level
              pop            X
              pop            Y
              call           next_level
              pop            X
              pop            Y
              call           next_level
              pop            X
              pop            Y
              shl            $2, N2
              jmp            coarse_level
Nlt64:        vmovapd        (X), er0
              vmovapd        32(X), er1
              vmovapd        64(X), er2
              vmovapd        96(X), er3
              vmovapd        (Y), ei0
              vmovapd        32(Y), ei1
              vmovapd        64(Y), ei2
              vmovapd        96(Y), ei3
              call           transpose_lo
              call           butterfly4
              call           transpose_lo
              vmovapd        TWO, two
              vmovapd        tw16_1r, cos1
              vmovapd        tw16_1i, sin1
              vmovapd        tw16_2r, cos2
              vmovapd        tw16_2i, sin2
              call           tw_butterfly4
              vmovapd        er0, (X)
              vmovapd        er3, 32(X)
              vmovapd        er1, 64(X)
              vmovapd        er2, 96(X)
              vmovapd        ei0, (Y)
              vmovapd        ei3, 32(Y)
              vmovapd        ei1, 64(Y)
              vmovapd        ei2, 96(Y)
              cmp            $8, N2
              jl             done
              vmovapd        128(X), er0
              vmovapd        160(X), er1
              vmovapd        192(X), er2
              vmovapd        224(X), er3
              vmovapd        128(Y), ei0
              vmovapd        160(Y), ei1
              vmovapd        192(Y), ei2
              vmovapd        224(Y), ei3
              call           transpose_lo
              call           butterfly4
              call           transpose_lo
              vmovapd        TWO, two
              vmovapd        tw16_1r, cos1
              vmovapd        tw16_1i, sin1
              vmovapd        tw16_2r, cos2
              vmovapd        tw16_2i, sin2
              call           tw_butterfly4
              vmovapd        er0, 128(X)
              vmovapd        er3, 160(X)
              vmovapd        er1, 192(X)
              vmovapd        er2, 224(X)
              vmovapd        ei0, 128(Y)
              vmovapd        ei3, 160(Y)
              vmovapd        ei1, 192(Y)
              vmovapd        ei2, 224(Y)
              vmovapd        (X), er0
              vmovapd        128(X), er2
              vmovapd        (Y), ei0
              vmovapd        128(Y), ei2
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
              vmovapd        er0, (X)
              vmovapd        er2, 128(X)
              vmovapd        ei0, (Y)
              vmovapd        ei2, 128(Y)
              vmovapd        32(X), er0
              vmovapd        160(X), er2
              vmovapd        32(Y), ei0
              vmovapd        160(Y), ei2
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
              vmovapd        er0, 32(X)
              vmovapd        er2, 160(X)
              vmovapd        ei0, 32(Y)
              vmovapd        ei2, 160(Y)
              vmovapd        64(X), er0
              vmovapd        192(X), er2
              vmovapd        64(Y), ei0
              vmovapd        192(Y), ei2
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
              vmovapd        er0, 64(X)
              vmovapd        er2, 192(X)
              vmovapd        ei0, 64(Y)
              vmovapd        ei2, 192(Y)
              vmovapd        96(X), er0
              vmovapd        224(X), er2
              vmovapd        96(Y), ei0
              vmovapd        224(Y), ei2
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
              vmovapd        er0, 96(X)
              vmovapd        er2, 224(X)
              vmovapd        ei0, 96(Y)
              vmovapd        ei2, 224(Y)
              jmp            done
coarse_level: bsf            N2, %rax
              add            $2, %rax
              shl            %rax
              mov            (Wptr, %rax, 8), Wr
              mov            8(Wptr, %rax, 8), Wi
              imul           $2, N2, %rcx
              imul           $3, N2, %rdx
              xor            k2, k2
k2_loop:      vmovapd        (Wr, k2, 8), cos1
              vmovapd        (Wi, k2, 8), sin1
              vmulpd         sin1, sin1, cos2
              vmulpd         cos1, sin1, sin2
              vfmsub231pd    cos1, cos1, cos2
              vfmadd231pd    sin1, cos1, sin2
              lea            (X, k2, 8), %rax
              lea            (Y, k2, 8), %rbx
              vmovapd        (%rax), er0
              vmovapd        (%rax, N2, 8), er2
              vmovapd        (%rax, %rcx, 8), er1
              vmovapd        (%rax, %rdx, 8), er3
              vmovapd        (%rbx), ei0
              vmovapd        (%rbx, N2, 8), ei2
              vmovapd        (%rbx, %rcx, 8), ei1
              vmovapd        (%rbx, %rdx, 8), ei3
              call           tw_butterfly4
              vmovapd        er0, (%rax)
              vmovapd        er3, (%rax, N2, 8)
              vmovapd        er1, (%rax, %rcx, 8)
              vmovapd        er2, (%rax, %rdx, 8)
              vmovapd        ei0, (%rbx)
              vmovapd        ei3, (%rbx, N2, 8)
              vmovapd        ei1, (%rbx, %rcx, 8)
              vmovapd        ei2, (%rbx, %rdx, 8)
              add            $4, k2
              cmp            k2, N2
              jne            k2_loop
done:         ret
              .align         32
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

