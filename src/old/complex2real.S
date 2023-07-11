#include      "ymm.h"
#define       X0             %r15
#define       Y0             %r14
#define       X1             %r13
#define       Y1             %r12
#define       N              %r11
#define       No4            %r10
#define       k2             %r9
#define       No4mk2         %r8
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
              push           %r10
              push           %r11
              call           twiddles_fwd_complex
              pop            %r11
              pop            %r10
              pop            %r8
              mov            %rax, Wr
              mov            %rdx, Wi
              vmovapd        TWO, two
              vmovq          (X0), ter0
              vmovq          (Y0), ter2
              vmovq          (X1), ter1
              vmovq          (Y1), ter3
              vaddsd         ter2, ter0, tei0
              vsubsd         ter2, ter0, tei2
              vaddsd         ter3, ter1, tei1
              vsubsd         ter1, ter3, tei3
              vaddsd         tei1, tei0, ter0
              vsubsd         tei1, tei0, ter2
              vmovq          ter0, (X0)
              vmovq          tei2, (Y0)
              vmovq          ter2, (X1)
              vmovq          tei3, (Y1)
              vmovq          (X0, No4, 4), tei0
              vmovq          (Y0, No4, 4), tei2
              vmovq          (X1, No4, 4), tei1
              vmovq          (Y1, No4, 4), tei3
              vaddsd         tei3, tei1, ter0
              vsubsd         tei3, tei1, ter2
              vmulsd         TW45, ter2, ter1
              vmulsd         TW45, ter0, ter3
              vmovq          tei0, ter0
              vmovq          tei2, ter2
              vaddsd         ter1, ter0, tei0
              vaddsd         ter3, ter2, tei3
              vsubsd         ter1, ter0, tei1
              vsubsd         ter3, ter2, tei2
              vmulsd         NONE, tei3, tei3
              vmovq          tei0, (X0, No4, 4)
              vmovq          tei1, (Y0, No4, 4)
              vmovq          tei2, (X1, No4, 4)
              vmovq          tei3, (Y1, No4, 4)
              mov            $1, k2
scalar_loop:  mov            No4, No4mk2
              sub            k2, No4mk2
              vmovq          (X0, k2, 8), ttr0
              vmovq          (X0, No4mk2, 8), ttr1
              vmovq          (Y0, k2, 8), tti0
              vmovq          (Y0, No4mk2, 8), tti1
              vmovq          (X1, k2, 8), ttr2
              vmovq          (X1, No4mk2, 8), ttr3
              vmovq          (Y1, k2, 8), tti2
              vmovq          (Y1, No4mk2, 8), tti3
              vaddsd         ttr1, ttr0, ter0
              vsubsd         tti1, tti0, tei0
              vsubsd         ttr0, ttr1, tei2
              vaddsd         tti1, tti0, ter2
              vaddsd         ttr3, ttr2, ter1
              vsubsd         tti3, tti2, tei1
              vsubsd         ttr2, ttr3, tei3
              vaddsd         tti3, tti2, ter3
              vmovq          (Wr, k2, 8), tcos1
              vmovq          (Wi, k2, 8), tsin1
              vmulpd         tsin1, tsin1, tcos2
              vmulpd         tcos1, tsin1, tsin2
              vfmsub231pd    tcos1, tcos1, tcos2
              vfmadd231pd    tsin1, tcos1, tsin2
              vmovq          TWO, ttwo
              call           tw_butterfly4_scalar
              vmulsd         NONE, tei1, tei1
              vmulsd         NONE, tei2, tei2
              vmovq          HALF, thalf
              vmulsd         thalf, ter0, ter0
              vmulsd         thalf, ter1, ter1
              vmulsd         thalf, tei0, tei0
              vmulsd         thalf, tei1, tei1
              vmulsd         thalf, ter2, ter2
              vmulsd         thalf, ter3, ter3
              vmulsd         thalf, tei2, tei2
              vmulsd         thalf, tei3, tei3
              vmovq          ter0, (X0, k2, 8)
              vmovq          ter3, (Y0, k2, 8)
              vmovq          tei1, (X1, k2, 8)
              vmovq          tei2, (Y1, k2, 8)
              vmovq          ter2, (X0, No4mk2, 8)
              vmovq          ter1, (Y0, No4mk2, 8)
              vmovq          tei3, (X1, No4mk2, 8)
              vmovq          tei0, (Y1, No4mk2, 8)
              inc            k2
              cmp            $4, k2
              jne            scalar_loop
              mov            $4, k2
simd_loop:    mov            No4, No4mk2
              sub            k2, No4mk2
              sub            $3, No4mk2
              vmovapd        (X0, k2, 8), tr0
              vpermpd        $27, (X0, No4mk2, 8), tr1
              vmovapd        (Y0, k2, 8), ti0
              vpermpd        $27, (Y0, No4mk2, 8), ti1
              vmovapd        (X1, k2, 8), tr2
              vpermpd        $27, (X1, No4mk2, 8), tr3
              vmovapd        (Y1, k2, 8), ti2
              vpermpd        $27, (Y1, No4mk2, 8), ti3
              vaddpd         tr1, tr0, er0
              vsubpd         ti1, ti0, ei0
              vsubpd         tr0, tr1, ei2
              vaddpd         ti1, ti0, er2
              vaddpd         tr3, tr2, er1
              vsubpd         ti3, ti2, ei1
              vsubpd         tr2, tr3, ei3
              vaddpd         ti3, ti2, er3
              vmovapd        (Wr, k2, 8), cos1
              vmovapd        (Wi, k2, 8), sin1
              vmulpd         sin1, sin1, cos2
              vmulpd         cos1, sin1, sin2
              vfmsub231pd    cos1, cos1, cos2
              vfmadd231pd    sin1, cos1, sin2
              vmovapd        TWO, two
              call           tw_butterfly4
              vpermpd        $27, ei0, ei0
              vpermpd        $27, ei3, ei3
              vpermpd        $27, er1, er1
              vpermpd        $27, er2, er2
              vmulpd         NONE, ei2, ei2
              vmulpd         NONE, ei1, ei1
              vmovapd        HALF, half
              vmulpd         half, er0, er0
              vmulpd         half, er1, er1
              vmulpd         half, ei0, ei0
              vmulpd         half, ei1, ei1
              vmulpd         half, er2, er2
              vmulpd         half, er3, er3
              vmulpd         half, ei2, ei2
              vmulpd         half, ei3, ei3
              vmovapd        er0, (X0, k2, 8)
              vmovapd        er3, (Y0, k2, 8)
              vmovapd        ei1, (X1, k2, 8)
              vmovapd        ei2, (Y1, k2, 8)
              vmovupd        er2, (X0, No4mk2, 8)
              vmovupd        er1, (Y0, No4mk2, 8)
              vmovupd        ei3, (X1, No4mk2, 8)
              vmovupd        ei0, (Y1, No4mk2, 8)
              add            $4, k2
              mov            No4, %rax
              sub            k2, %rax
              cmp            k2, %rax
              jg             simd_loop
done:         pop            %r15
              pop            %r14
              pop            %r13
              pop            %r12
              pop            %rbx
              ret
              .align         8
TW45:         .double        0.70710678118654752440
              .align         32
TWO:          .double        2.0
              .double        2.0
              .double        2.0
              .double        2.0
HALF:         .double        0.5
              .double        0.5
              .double        0.5
              .double        0.5
NONE:         .double        -1.0
              .double        -1.0
              .double        -1.0
              .double        -1.0
