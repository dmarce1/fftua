#include      "ymm.h"
#define       X0             %r15
#define       Y0             %r14
#define       X1             %r13
#define       Y1             %r12
#define       N              %r11
#define       No4            %r10
#define       k2             %r9
#define       Nmk2           %r8
#define       Wr             %rsi
#define       Wi             %rdi

              .global        complex2real


              .text

complex2real: push           %rbx
              push           %r12
              push           %r13
              push           %r14
              push           %r15
              mov            %rdi, X0
              mov            %rsi, Y0
              mov            %rdx, X1
              mov            %rcx, Y1
              mov            %r8, N
              mov            N, No4
              shr            $2, No4
              mov            N, %rdi
              push           %r8
              call           twiddles_fwd_complex
              pop            %r8
              mov            %rax, Wr
              mov            %rdx, Wi
              vmovapd        HALF, half
              mov            $4, k2
simd_loop:    mov            N, Nmk2
              sub            k2, Nmk2
              sub            $3, Nmk2
              vmovapd        (Wr, k2, 8), cos1
              vmovapd        (Wi, k2, 8), sin1
              vmulpd         sin1, sin1, cos2
              vmulpd         cos1, sin1, sin2
              vfmsub231pd    cos1, cos1, cos2
              vfmadd231pd    sin1, cos1, sin2
              vmovapd        (X0, k2, 8), er0
              vpermpd        $27, (X0, Nmk2, 8), er1
              vmovapd        (Y0, k2, 8), ei0
              vpermpd        $27, (Y0, Nmk2, 8), ei1
              vmovapd        (X1, k2, 8), er2
              vpermpd        $27, (X1, Nmk2, 8), er3
              vmovapd        (Y1, k2, 8), ei2
              vpermpd        $27, (Y1, Nmk2, 8), ei3
              vaddpd         er1, er0, tr0
              vsubpd         ei1, ei0, ti0
              vsubpd         er0, er1, ti1
              vaddpd         ei1, ei0, tr1
              vmulpd         half, tr0, er0
              vmulpd         half, tr1, er1
              vmulpd         half, ti0, ei0
              vmulpd         half, ti1, ei1
              call           tw_butterfly4
              vpermpd        $27, ei0, ei0
              vpermpd        $27, ei1, ei1
              vpermpd        $27, ei2, ei2
              vpermpd        $27, ei3, ei3
              vmulpd         NONE, ei3, ei3
              vmulpd         NONE, ei2, ei2
              vmovapd        er0, (X0, k2, 8)
              vmovapd        er3, (X0, Nmk2, 8)
              vmovapd        er1, (Y0, k2, 8)
              vmovapd        er2, (Y0, Nmk2, 8)
              vmovapd        ei2, (X1, k2, 8)
              vmovapd        ei1, (X1, Nmk2, 8)
              vmovapd        ei3, (Y1, k2, 8)
              vmovapd        ei0, (Y1, Nmk2, 8)
              add            $4, k2
              cmp            k2, No4
              jne            simd_loop
              pop            %r15
              pop            %r14
              pop            %r13
              pop            %r12
              pop            %rbx
              ret
              .align         32
HALF:         .double        0.5
              .double        0.5
              .double        0.5
              .double        0.5
NONE:         .double        -1.0
              .double        -1.0
              .double        -1.0
              .double        -1.0
