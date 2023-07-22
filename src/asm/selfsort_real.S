#define       X              %r15
#define       N              %r14
#define       Wr             %r13
#define       Wi             %r12
#define       I              %r11
#define       M              %r10
#define       k2lo           %r9
#define       k2hi           %r8
#define       er0            %ymm0
#define       ei0            %ymm1
#define       er1            %ymm2
#define       ei1            %ymm3
#define       er2            %ymm4
#define       ei2            %ymm5
#define       er3            %ymm6
#define       ei3            %ymm7
#define       cos1           %ymm8
#define       sin1           %ymm9
#define       cos2           %ymm10
#define       sin2           %ymm11
#define       tr             %ymm12
#define       ti             %ymm13
#define       two            %ymm14
#define       ter0           %xmm0
#define       tei0           %xmm1
#define       ter1           %xmm2
#define       tei1           %xmm3
#define       ter2           %xmm4
#define       tei2           %xmm5
#define       ter3           %xmm6
#define       tei3           %xmm7
#define       tcos1          %xmm8
#define       tsin1          %xmm9
#define       tcos2          %xmm10
#define       tsin2          %xmm11
#define       ttr            %xmm12
#define       tti            %xmm13
#define       ttwo           %xmm14
#define       SIMD_SIZE      $4

              .global        fft_selfsort_real2


              .text

fft_selfsort_real2:
              push           %rbx
              push           %r15
              push           %r14
              push           %r13
              push           %r12
              mov            %rdi, X
              mov            %rsi, N
              mov            N, %rdi
              call           get_forward_twiddles
              mov            %rax, Wr
              mov            %rdx, Wi
              call           trivial_butterflies
              call           transpose_outer
              call           final_butterflies
              pop            %r12
              pop            %r13
              pop            %r14
              pop            %r15
              pop            %rbx
              ret

next_level:


trivial_butterflies:
              mov            N, M
              shr            $2, M
              xor            I, I
L1000:        lea            (X, I, 8), %rax
              lea            (%rax, M, 8), %rbx
              lea            (%rbx, M, 8), %rcx
              lea            (%rcx, M, 8), %rdx
              vmovapd        (%rax), er0
              vmovapd        (%rbx), er1
              vmovapd        (%rcx), er2
              vmovapd        (%rdx), er3
              call           butterfly04_simd
              vmovapd        er0, (%rax)
              vmovapd        ei2, (%rbx)
              vmovapd        ei3, (%rcx)
              vmovapd        er2, (%rdx)
              add            SIMD_SIZE, I
              cmp            I, M
              jne            L1000
              ret

transpose_outer:
              mov            N, M
              shr            $2, M
              xor            I, I
L2000:        mov            M, %rdi
              sub            I, %rdi
              sub            I, %rdi
              sub            $3, %rdi
              lea            (X, I, 8), %rax
              lea            (%rax, M, 8), %rbx
              lea            (%rbx, M, 8), %rcx
              lea            (%rcx, M, 8), %rdx
              vmovapd        (%rax), er0
              vmovapd        (%rbx), er1
              vmovapd        (%rcx), er2
              vmovapd        (%rdx), er3
              vunpcklpd      er1, er0, %ymm12
              vunpckhpd      er1, er0, %ymm13
              vunpcklpd      er2, er3, %ymm14
              vunpckhpd      er2, er3, %ymm15
              vperm2f128     $0x20, %ymm14, %ymm12, er0
              vperm2f128     $0x20, %ymm15, %ymm13, er1
              vperm2f128     $0x31, %ymm14, %ymm12, er2
              vperm2f128     $0x31, %ymm15, %ymm13, er3
              vmovapd        er0, (%rax)
              vmovapd        er1, (%rbx)
              vmovapd        er2, (%rcx)
              vmovapd        er3, (%rdx)
              add            SIMD_SIZE, I
              cmp            I, M
              jne            L2000
              ret

final_butterflies:
              mov            N, M
              shr            $2, M
              lea            (X), %rax
              lea            (%rax, M, 8), %rbx
              lea            (%rbx, M, 8), %rcx
              lea            (%rcx, M, 8), %rdx
              vmovq          (%rax), ter0
              vmovq          (%rbx), ter1
              vmovq          (%rcx), ter2
              vmovq          (%rdx), ter3
              call           butterfly04_scalar
              vmovq          ter0, (%rax)
              vmovq          tei2, (%rbx)
              vmovq          ter2, (%rcx)
              vmovq          tei3, (%rdx)
              lea            (X, M, 4), %rax
              lea            (%rax, M, 8), %rbx
              lea            (%rbx, M, 8), %rcx
              lea            (%rcx, M, 8), %rdx
              vmovq          (%rax), ter0
              vmovq          (%rbx), ter1
              vmovq          (%rcx), ter2
              vmovq          (%rdx), ter3
              call           butterflyNy4_scalar
              vmovq          ter0, (%rax)
              vmovq          ter1, (%rbx)
              vmovq          ter2, (%rcx)
              vmovq          ter3, (%rdx)
              mov            M, %rsi
              shr            %rsi
              mov            $4, %rax
              cmp            $4, %rsi
              cmovg          %rax, %rsi
              cmp            $1, %rsi
              jle            final_done
              mov            $1, I
L3500:        vmovq          (Wr, I, 8), tcos1
              vmovq          (Wi, I, 8), tsin1
              vmulsd         tsin1, tsin1, tcos2
              vmulsd         tcos1, tsin1, tsin2
              vfmsub231sd    tcos1, tcos1, tcos2
              vfmadd231sd    tsin1, tcos1, tsin2
              mov            M, %rdi
              sub            I, %rdi
              sub            I, %rdi
              lea            (X, I, 8), %rax
              lea            (%rax, M, 8), %rbx
              lea            (%rbx, M, 8), %rcx
              lea            (%rcx, M, 8), %rdx
              vmovq          (%rax), ter0
              vmovq          (%rbx), ter1
              vmovq          (%rcx), ter2
              vmovq          (%rdx), ter3
              vmovq          (%rax, %rdi, 8), tei0
              vmovq          (%rbx, %rdi, 8), tei1
              vmovq          (%rcx, %rdi, 8), tei2
              vmovq          (%rdx, %rdi, 8), tei3
              call           butterfly4_scalar
              vmovq          ter0, (%rax)
              vmovq          ter2, (%rax, %rdi, 8)
              vmovq          ter3, (%rbx)
              vmovq          ter1, (%rbx, %rdi, 8)
              vmovq          tei1, (%rcx)
              vmovq          tei3, (%rcx, %rdi, 8)
              vmovq          tei2, (%rdx)
              vmovq          tei0, (%rdx, %rdi, 8)
              inc            I
              cmp            I, %rsi
              jne            L3500
              mov            M, %rsi
              shr            %rsi
              mov            $4, I
L4000:        cmp            I, %rsi
              jle            final_done
              vmovapd        (Wr, I, 8), cos1
              vmovapd        (Wi, I, 8), sin1
              vmulpd         sin1, sin1, cos2
              vmulpd         cos1, sin1, sin2
              vfmsub231pd    cos1, cos1, cos2
              vfmadd231pd    sin1, cos1, sin2
              mov            M, %rdi
              sub            I, %rdi
              sub            I, %rdi
              sub            $3, %rdi
              lea            (X, I, 8), %rax
              lea            (%rax, M, 8), %rbx
              lea            (%rbx, M, 8), %rcx
              lea            (%rcx, M, 8), %rdx
              vmovapd        (%rax), er0
              vmovapd        (%rbx), er1
              vmovapd        (%rcx), er2
              vmovapd        (%rdx), er3
              vpermpd        $27, (%rax, %rdi, 8), ei0
              vpermpd        $27, (%rbx, %rdi, 8), ei1
              vpermpd        $27, (%rcx, %rdi, 8), ei2
              vpermpd        $27, (%rdx, %rdi, 8), ei3
              call           butterfly4_simd
              vpermpd        $27, er2, er2
              vpermpd        $27, er1, er1
              vpermpd        $27, ei0, ei0
              vpermpd        $27, ei3, ei3
              vmovapd        er0, (%rax)
              vmovapd        er2, (%rax, %rdi, 8)
              vmovapd        er3, (%rbx)
              vmovapd        er1, (%rbx, %rdi, 8)
              vmovapd        ei1, (%rcx)
              vmovapd        ei3, (%rcx, %rdi, 8)
              vmovapd        ei2, (%rdx)
              vmovapd        ei0, (%rdx, %rdi, 8)
              add            SIMD_SIZE, I
              jmp            L4000
final_done:   ret

butterfly04_simd:
              vaddpd         er2, er0, ei0
              vsubpd         er2, er0, ei2
              vaddpd         er3, er1, ei1
              vsubpd         er1, er3, ei3
              vaddpd         ei1, ei0, er0
              vsubpd         ei1, ei0, er2
              ret

butterfly04_scalar:
              vaddpd         ter2, ter0, tei0
              vsubpd         ter2, ter0, tei2
              vaddpd         ter3, ter1, tei1
              vsubpd         ter1, ter3, tei3
              vaddpd         tei1, tei0, ter0
              vsubpd         tei1, tei0, ter2
              ret

butterflyNy4_simd:
              vaddpd         er3, er1, ei0
              vsubpd         er3, er1, ei2
              vmulpd         TW45, ei2, ei1
              vmulpd         TW45, ei0, ei3
              vmovapd        er0, ei0
              vmovapd        er2, ei2
              vaddpd         ei1, ei0, er0
              vaddpd         ei3, ei2, er3
              vsubpd         ei1, ei0, er1
              vsubpd         ei3, ei2, er2
              vmulpd         NONE, er3, er3
              ret

butterflyNy4_scalar:
              vaddsd         ter3, ter1, tei0
              vsubsd         ter3, ter1, tei2
              vmulsd         TW45, tei2, tei1
              vmulsd         TW45, tei0, tei3
              vmovq          ter0, tei0
              vmovq          ter2, tei2
              vaddsd         tei1, tei0, ter0
              vaddsd         tei3, tei2, ter3
              vsubsd         tei1, tei0, ter1
              vsubsd         tei3, tei2, ter2
              vmulsd         NONE, ter3, ter3
              ret

butterfly4_simd:
              vmovapd        TWO, two
              vmovapd        er0, tr
              vmovapd        ei0, ti
              vfmadd231pd    sin2, ei2, tr
              vfnmadd231pd   sin2, er2, ti
              vfnmadd132pd   cos2, tr, er2
              vfnmadd132pd   cos2, ti, ei2
              vfmsub132pd    two, er2, er0
              vfmsub132pd    two, ei2, ei0
              vmovapd        er1, tr
              vmovapd        ei1, ti
              vfmadd231pd    sin2, ei3, tr
              vfnmadd231pd   sin2, er3, ti
              vfnmadd132pd   cos2, tr, er3
              vfnmadd132pd   cos2, ti, ei3
              vfmsub132pd    two, er3, er1
              vfmsub132pd    two, ei3, ei1
              vmovapd        er0, tr
              vmovapd        ei0, ti
              vfmadd231pd    sin1, ei1, tr
              vfnmadd231pd   sin1, er1, ti
              vfnmadd132pd   cos1, tr, er1
              vfnmadd132pd   cos1, ti, ei1
              vfmsub132pd    two, er1, er0
              vfmsub132pd    two, ei1, ei0
              vmovapd        er2, tr
              vmovapd        ei2, ti
              vfmadd231pd    cos1, ei3, tr
              vfnmadd231pd   cos1, er3, ti
              vfmadd132pd    sin1, tr, er3
              vfmadd132pd    sin1, ti, ei3
              vfmsub132pd    two, er3, er2
              vfmsub132pd    two, ei3, ei2
              vmulpd         NONE, ei1, ei1
              vmulpd         NONE, ei2, ei2
              ret

butterfly4_scalar:
              vmovq          TWO, ttwo
              vmovq          ter0, ttr
              vmovq          tei0, tti
              vfmadd231sd    tsin2, tei2, ttr
              vfnmadd231sd   tsin2, ter2, tti
              vfnmadd132sd   tcos2, ttr, ter2
              vfnmadd132sd   tcos2, tti, tei2
              vfmsub132sd    ttwo, ter2, ter0
              vfmsub132sd    ttwo, tei2, tei0
              vmovq          ter1, ttr
              vmovq          tei1, tti
              vfmadd231sd    tsin2, tei3, ttr
              vfnmadd231sd   tsin2, ter3, tti
              vfnmadd132sd   tcos2, ttr, ter3
              vfnmadd132sd   tcos2, tti, tei3
              vfmsub132sd    ttwo, ter3, ter1
              vfmsub132sd    ttwo, tei3, tei1
              vmovq          ter0, ttr
              vmovq          tei0, tti
              vfmadd231sd    tsin1, tei1, ttr
              vfnmadd231sd   tsin1, ter1, tti
              vfnmadd132sd   tcos1, ttr, ter1
              vfnmadd132sd   tcos1, tti, tei1
              vfmsub132sd    ttwo, ter1, ter0
              vfmsub132sd    ttwo, tei1, tei0
              vmovq          ter2, ttr
              vmovq          tei2, tti
              vfmadd231sd    tcos1, tei3, ttr
              vfnmadd231sd   tcos1, ter3, tti
              vfmadd132sd    tsin1, ttr, ter3
              vfmadd132sd    tsin1, tti, tei3
              vfmsub132sd    ttwo, ter3, ter2
              vfmsub132sd    ttwo, tei3, tei2
              vmulsd         NONE, tei1, tei1
              vmulsd         NONE, tei2, tei2
              ret

              .align         32
TWO:          .double        2.0
              .double        2.0
              .double        2.0
              .double        2.0
TW45:         .double        0.70710678118654
              .double        0.70710678118654
              .double        0.70710678118654
              .double        0.70710678118654
NONE:         .double        -1.0
              .double        -1.0
              .double        -1.0
              .double        -1.0

