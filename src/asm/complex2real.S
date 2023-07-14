#define       er0            %ymm0
#define       ei0            %ymm1
#define       er1            %ymm2
#define       ei1            %ymm3
#define       er2            %ymm4
#define       ei2            %ymm5
#define       er3            %ymm6
#define       ei3            %ymm7
#define       tr0            %ymm8
#define       tr1            %ymm9
#define       tr2            %ymm10
#define       tr3            %ymm11
#define       ti0            %ymm12
#define       ti1            %ymm13
#define       ti2            %ymm14
#define       ti3            %ymm15
#define       cos1           %ymm8
#define       sin1           %ymm9
#define       cos2           %ymm10
#define       sin2           %ymm11
#define       tr             %ymm12
#define       ti             %ymm13
#define       two            %ymm14
#define       half           %ymm15
#define       ter0           %xmm0
#define       tei0           %xmm1
#define       ter1           %xmm2
#define       tei1           %xmm3
#define       ter2           %xmm4
#define       tei2           %xmm5
#define       ter3           %xmm6
#define       tei3           %xmm7
#define       ttr0           %xmm8
#define       ttr1           %xmm9
#define       ttr2           %xmm10
#define       ttr3           %xmm11
#define       tti0           %xmm12
#define       tti1           %xmm13
#define       tti2           %xmm14
#define       tti3           %xmm15
#define       tcos1          %xmm8
#define       tsin1          %xmm9
#define       tcos2          %xmm10
#define       tsin2          %xmm11
#define       ttr            %xmm12
#define       tti            %xmm13
#define       ttwo           %xmm14
#define       thalf          %xmm15
#define       X              %r15
#define       Wr             %r14
#define       Wi             %r13
#define       N              %r12
#define       k2             %r11
#define       M              %r10
#define       I              %r9


              .global        complex2real


              .text

complex2real: push           %rbx
              push           %r15
              push           %r13
              push           %r14
              push           %r12
              mov            %rdi, X
              mov            %rsi, N
              mov            N, %rdi
              call           get_forward_twiddles
              mov            %rax, Wr
              mov            %rdx, Wi
              mov            N, M
              shr            M
              lea            (X), %rax
              lea            (%rax, M, 4), %rbx
              lea            (%rbx, M, 4), %rcx
              lea            (%rcx, M, 4), %rdx
              vmovq          (%rax), ter0
              vmovq          (%rbx), ter1
              vmovq          (%rcx), tei0
              vmovq          (%rdx), tei1
              vaddpd         tei0, ter0, ttr0
              vsubpd         tei0, ter0, ttr1
              vmovq          ter1, ttr
              vmulsd         NONE, tei1, tti
              vmovq          ttr0, (%rax)
              vmovq          ttr, (%rbx)
              vmovq          ttr1, (%rcx)
              vmovq          tti, (%rdx)
              mov            $1, k2
final_scalar: mov            M, %rdi
              sub            k2, %rdi
              sub            k2, %rdi
              lea            (X, k2, 8), %rax
              lea            (%rax, M, 8), %rdx
              vmovq          (%rax), ter0
              vmovq          (%rax, %rdi, 8), ter1
              vmovq          (%rdx), tei0
              vmovq          (%rdx, %rdi, 8), tei1
              vaddsd         ter1, ter0, ttr0
              vsubsd         tei1, tei0, tti0
              vaddsd         tei0, tei1, ttr1
              vsubsd         ter0, ter1, tti1
              vmovq          HALF, thalf
              vmulsd         thalf, ttr0, ter0
              vmulsd         thalf, ttr1, ter1
              vmulsd         thalf, tti0, tei0
              vmulsd         thalf, tti1, tei1
              vmovq          (Wr, k2, 8), tcos1
              vmovq          (Wi, k2, 8), tsin1
              vmovq          TWO, ttwo
              vmovq          ter0, ttr
              vmovq          tei0, tti
              vfmadd231sd    tsin1, tei1, ttr
              vfnmadd231sd   tsin1, ter1, tti
              vfnmadd132sd   tcos1, ttr, ter1
              vfnmadd132sd   tcos1, tti, tei1
              vfmsub132sd    ttwo, ter1, ter0
              vfmsub132sd    ttwo, tei1, tei0
              vmulsd         NONE, tei1, tei1
              vmovq          ter0, (%rax)
              vmovq          ter1, (%rax, %rdi, 8)
              vmovq          tei1, (%rdx)
              vmovq          tei0, (%rdx, %rdi, 8)
              inc            k2
              cmp            $4, k2
              jne            final_scalar
              mov            $4, k2
final_simd:   mov            M, %rdi
              sub            k2, %rdi
              sub            k2, %rdi
              sub            $3, %rdi
              lea            (X, k2, 8), %rax
              lea            (%rax, M, 8), %rdx
              vmovapd        (%rax), er0
              vmovapd        (%rdx), ei0
              vpermpd        $27, (%rax, %rdi, 8), er1
              vpermpd        $27, (%rdx, %rdi, 8), ei1
              vaddpd         er1, er0, tr0
              vsubpd         ei1, ei0, ti0
              vaddpd         ei0, ei1, tr1
              vsubpd         er0, er1, ti1
              vmovapd        HALF, half
              vmulpd         half, tr0, er0
              vmulpd         half, tr1, er1
              vmulpd         half, ti0, ei0
              vmulpd         half, ti1, ei1
              vmovapd        (Wr, k2, 8), cos1
              vmovapd        (Wi, k2, 8), sin1
              vmovapd        TWO, two
              vmovapd        er0, tr
              vmovapd        ei0, ti
              vfmadd231pd    sin1, ei1, tr
              vfnmadd231pd   sin1, er1, ti
              vfnmadd132pd   cos1, tr, er1
              vfnmadd132pd   cos1, ti, ei1
              vfmsub132pd    two, er1, er0
              vfmsub132pd    two, ei1, ei0
              vmulpd         NONE, ei1, ei1
              vpermpd        $27, er1, er1
              vpermpd        $27, ei0, ei0
              vmovapd        er0, (%rax)
              vmovapd        ei1, (%rdx)
              vmovupd        er1, (%rax, %rdi, 8)
              vmovupd        ei0, (%rdx, %rdi, 8)
              add            $4, k2
              mov            M, %rax
              shr            %rax
              cmp            k2, %rax
              jne            final_simd
              pop            %r12
              pop            %r14
              pop            %r13
              pop            %r15
              pop            %rbx
              ret
