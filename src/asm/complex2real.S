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
#define       half           %ymm14
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
#define       thalf          %xmm14
#define       X              %r15
#define       N              %r14
#define       I              %r13
#define       N2             %r12
#define       k2             %r11
#define       Wr             %r10
#define       Wi             %r9
#define       N2o2           %r8

//ur0 ur1 ur2 ur3 ui0 ui1 ui2 ui3
// e0  o0  e1  o1  e2  o2  e3  o3
              .global        complex2real

              .text

complex2real: push           %r15
              push           %r14
              push           %r13
              push           %r12
              push           %rbx
              mov            %rdi, X
              mov            %rsi, N
              mov            N, N2
              shr            $3, N2
              mov            N, %rdi
              call           get_forward_twiddles
              push           %rdx
              push           %rax
              mov            N, %rdi
              shr            %rdi
              call           get_forward_twiddles
              mov            %rax, Wr
              mov            %rdx, Wi
              mov            N2, %rdi
              shl            $2, %rdi
              xor            k2, k2
              vmovapd        TWO, two
radix4:       vmovapd        (Wr, k2, 8), cos1
              vmovapd        (Wi, k2, 8), sin1
              vmulpd         sin1, sin1, cos2
              vmulpd         cos1, sin1, sin2
              vfmsub231pd    cos1, cos1, cos2
              vfmadd231pd    sin1, cos1, sin2
              lea            (X, k2, 8), %rax
              lea            (%rax, N2, 8), %rbx
              lea            (%rbx, N2, 8), %rcx
              lea            (%rcx, N2, 8), %rdx
              vmovapd        (%rax), er0
              vmovapd        (%rbx), ei0
              vmovapd        (%rcx), er2
              vmovapd        (%rdx), ei2
              vmovapd        (%rax, %rdi, 8), er1
              vmovapd        (%rbx, %rdi, 8), ei1
              vmovapd        (%rcx, %rdi, 8), er3
              vmovapd        (%rdx, %rdi, 8), ei3
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
              vmovapd        er0, (%rax)
              vmovapd        er3, (%rbx)
              vmovapd        er1, (%rcx)
              vmovapd        er2, (%rdx)
              vmovapd        ei0, (%rax, %rdi, 8)
              vmovapd        ei3, (%rbx, %rdi, 8)
              vmovapd        ei1, (%rcx, %rdi, 8)
              vmovapd        ei2, (%rdx, %rdi, 8)
              add            $4, k2
              cmp            k2, N2
              jne            radix4
              pop            Wr
              pop            Wi
              mov            N, N2
              shr            N2
              lea            (X), %rax
              lea            (%rax, N2, 4), %rbx
              lea            (%rbx, N2, 4), %rcx
              lea            (%rcx, N2, 4), %rdx
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
radix2_sclar: mov            N2, %rdi
              sub            k2, %rdi
              sub            k2, %rdi
              lea            (X, k2, 8), %rax
              lea            (%rax, N2, 8), %rcx
              lea            (%rax, %rdi, 8), %rbx
              lea            (%rcx, %rdi, 8), %rdx
              vmovq          (%rax), ter0
              vmovq          (%rbx), ter1
              vmovq          (%rcx), tei0
              vmovq          (%rdx), tei1
              vaddsd         ter1, ter0, ttr0
              vsubsd         tei1, tei0, tti0
              vsubsd         ter0, ter1, tti1
              vaddsd         tei1, tei0, ttr1
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
              vmovq          ter1, (%rbx)
              vmovq          tei1, (%rcx)
              vmovq          tei0, (%rdx)
              inc            k2
              cmp            $4, k2
              jne            radix2_sclar
              mov            $4, k2
              mov            N2, N2o2
              shr            N2o2
radix2:       mov            N2, %rdi
              sub            k2, %rdi
              sub            k2, %rdi
              sub            $3, %rdi
              lea            (X, k2, 8), %rax
              lea            (%rax, %rdi, 8), %rbx
              lea            (%rax, N2, 8), %rcx
              lea            (%rcx, %rdi, 8), %rdx
              vmovapd        (%rax), er0
              vpermpd        $27, (%rbx), er1
              vmovapd        (%rcx), ei0
              vpermpd        $27, (%rdx), ei1
              vaddpd         er1, er0, tr0
              vsubpd         ei1, ei0, ti0
              vaddpd         ei1, ei0, tr1
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
              vmovupd        er1, (%rbx)
              vmovapd        ei1, (%rcx)
              vmovupd        ei0, (%rdx)
              add            $4, k2
              cmp            k2, N2o2
              jne            radix2
done:         pop            %rbx
              pop            %r12
              pop            %r13
              pop            %r14
              pop            %r15
              ret
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
ZERO:         .double         0.0
              .double         0.0
              .double         0.0
              .double         0.0
ONE:          .double         1.0
              .double         1.0
              .double         1.0
              .double         1.0
