#define N1            $4
#define SIMD_SIZE     $4
#define X             %r8
#define Y             %r9
#define W             %r10
#define Wv            %r11
#define NHI           %r12
#define N             %r13
#define N2            %r14
#define k2            %r15
#define NLO           %r15
#define I             %rax
#define J             %rbx
#define I4            %rcx
#define J4            %rdx
#define NMID          %r11
#define K             %rdi
#define L             %r10
#define M             %rsi
#define tmp           %ymm0
#define cos1          %ymm0
#define sin1          %ymm1
#define cos2          %ymm2
#define sin2          %ymm3
#define cos3          %ymm4
#define sin3          %ymm5
#define tr            %ymm6
#define ti            %ymm7
#define tr0           %ymm0
#define tr1           %ymm1
#define tr2           %ymm2
#define tr3           %ymm3
#define ti0           %ymm4
#define ti1           %ymm5
#define ti2           %ymm6
#define ti3           %ymm7
#define tr4           %ymm4
#define tr5           %ymm5
#define tr6           %ymm6
#define tr7           %ymm7
#define er0           %ymm8
#define er1           %ymm9
#define er2           %ymm10
#define er3           %ymm11
#define ei0           %ymm12
#define ei1           %ymm13
#define ei2           %ymm14
#define ei3           %ymm15
#define er4           %ymm12
#define er5           %ymm13
#define er6           %ymm14
#define er7           %ymm15
#define xcos1         %xmm0
#define xsin1         %xmm1
#define xcos2         %xmm2
#define xsin2         %xmm3
#define xcos3         %xmm4
#define xsin3         %xmm5
#define xtr           %xmm6
#define xti           %xmm7
#define xtr0          %xmm0
#define xtr1          %xmm1
#define xtr2          %xmm2
#define xtr3          %xmm3
#define xti0          %xmm4
#define xti1          %xmm5
#define xti2          %xmm6
#define xti3          %xmm7
#define xtr4          %xmm4
#define xtr5          %xmm5
#define xtr6          %xmm6
#define xtr7          %xmm7
#define xer0          %xmm8
#define xer1          %xmm9
#define xer2          %xmm10
#define xer3          %xmm11
#define xei0          %xmm12
#define xei1          %xmm13
#define xei2          %xmm14
#define xei3          %xmm15
#define xer4          %xmm12
#define xer5          %xmm13
#define xer6          %xmm14
#define xer7          %xmm15
#define pmask         $27

       .global        fft_recursive

       .data
       .align         32
iota:
       .quad          0
       .quad          2
       .quad          4
       .quad          6
ones:
       .quad          0xffffffffffffffff
       .quad          0xffffffffffffffff
       .quad          0xffffffffffffffff
       .quad          0xffffffffffffffff
C:
C0:    .double        0.0
C1:    .double        0.0
C2:    .double        0.0
C3:    .double        0.0
nC:
nC0:   .double        0.0
nC1:   .double        0.0
nC2:   .double        0.0
nC3:   .double        0.0
tmpr:
tmpr0:  .double       0.0
tmpr1:  .double       0.0
tmpr2:  .double       0.0
tmpr3:  .double       0.0
tmpi:
tmpi0:  .double       0.0
tmpi1:  .double       0.0
tmpi2:  .double       0.0
tmpi3:  .double       0.0
Z:
Z0:    .double        0.0
Z1:    .double        0.0
Z2:    .double        0.0
Z3:    .double        0.0
two:
       .double        2.0
       .double        2.0
       .double        2.0
       .double        2.0


       .text

fft_recursive:
       push           %rbx
       push           %r12
       push           %r13
       push           %r14
       push           %r15
       mov            %r8, N
       mov            %rdi, Y
       mov            %rsi, X
       mov            %rdx, W
       mov            %rcx, Wv
       shr            $2, N
       emms
       fld1
       fst            C0
       fst            C1
       fst            C2
       fst            C3
       vmovapd        C, tr
       vaddpd         tr, tr, ti
       vdivpd         ti, tr, tr
       vsqrtpd        tr, tr
       vmovapd        tr, C
       vmovapd        Z, ti
       vsubpd         C, ti, tr
       vmovapd        tr, nC
       call           scramble_hi
       call           transpose
       call           scramble_hi
       mov            $1, NHI
       mov            N, N2
       shr            $2, N2
       call           recurse
       shl            $2, N
       mov            N, N2
       shr            $2, N2
       call           finish
       pop            %r15
       pop            %r14
       pop            %r13
       pop            %r12
       pop            %rbx
       emms
       ret

recurse:
       cmp            $8, N
       jg             rcalls
       jl             skip_rcalls
       mov            $32, %rsi
       lea            (X), %rax
       lea            (%rax, %rsi), %rbx
       lea            (%rbx, %rsi), %rcx
       lea            (%rcx, %rsi), %rdx
       vmovapd        (%rax), er0
       vmovapd        (%rbx), er1
       vmovapd        (%rcx), er2
       vmovapd        (%rdx), er3
       vaddpd         er1, er0, tr0
       vsubpd         er1, er0, tr1
       vaddpd         er3, er2, tr2
       vsubpd         er3, er2, tr3
       vmovapd        tr0, (%rax)
       vmovapd        tr1, (%rbx)
       vmovapd        tr2, (%rcx)
       vmovapd        tr3, (%rdx)
       lea            (%rdx, %rsi), %rax
       lea            (%rax, %rsi), %rbx
       lea            (%rbx, %rsi), %rcx
       lea            (%rcx, %rsi), %rdx
       vmovapd        (%rax), er0
       vmovapd        (%rbx), er1
       vmovapd        (%rcx), er2
       vmovapd        (%rdx), er3
       vaddpd         er1, er0, tr0
       vsubpd         er1, er0, tr1
       vaddpd         er3, er2, tr2
       vsubpd         er3, er2, tr3
       vmovapd        tr0, (%rax)
       vmovapd        tr1, (%rbx)
       vmovapd        tr2, (%rcx)
       vmovapd        tr3, (%rdx)
       jmp            skip_rcalls
rcalls:
       push           X
       lea            (X, N, 8), %rax
       shr            $2, N2
       lea            (%rax, N, 8), %rbx
       shl            $2, NHI
       push           %rax
       lea            (%rbx, N, 8), X
       shr            $2, N
       push           %rbx
       call           recurse
       pop            X
       call           recurse
       pop            X
       call           recurse
       pop            X
       call           recurse
       shr            $2, NHI
       shl            $2, N
       shl            $2, N2
skip_rcalls:
       lea            (X), %rax
       lea            (%rax, N, 8), %rbx
       lea            (%rbx, N, 8), %rcx
       lea            (%rcx, N, 8), %rdx
       vmovapd        (%rax), er0
       vmovapd        (%rbx), er2
       vmovapd        (%rcx), er1
       vmovapd        (%rdx), er3
       vaddpd         er2, er0, tr0
       vsubpd         er2, er0, tr1
       vaddpd         er3, er1, tr2
       vsubpd         er1, er3, tr3
       vaddpd         tr2, tr0, er0
       vsubpd         tr2, tr0, er1
       vmovapd        er0, (%rax)
       vmovapd        tr1, (%rbx)
       vmovapd        er1, (%rcx)
       vmovapd        tr3, (%rdx)
       cmp            $1, N2
       jle            k2neqN2o2
       shl            N2
       lea            (X, N2, 8), %rax
       shr            N2
       lea            (%rax, N, 8), %rbx
       lea            (%rbx, N, 8), %rcx
       lea            (%rcx, N, 8), %rdx
       vmovapd        (%rax), tr0
       vmovapd        (%rbx), tr3
       vmovapd        (%rcx), tr1
       vmovapd        (%rdx), tr2
       vsubpd         tr2, tr1, er1
       vaddpd         tr2, tr1, er2
       vmulpd         C, er1, tr1
       vmulpd         nC, er2, tr2
       vaddpd         tr1, tr0, er0
       vsubpd         tr3, tr2, er1
       vsubpd         tr1, tr0, er2
       vaddpd         tr2, tr3, er3
       vmovapd        er1, (%rdx)
       vmovapd        er3, (%rcx)
       vmovapd        er2, (%rbx)
       vmovapd        er0, (%rax)
k2neqN2o2:
       xor            k2, k2
k2loopb:
       mov            N2, %rcx
       inc            k2
       shr            %rcx
       cmp            k2, %rcx
       jle            k2loope
       mov            NHI, %rax
       mul            k2
       shl            $3, %rax
       mov            %rax, %rcx
       mov            %rax, %rbx
       add            %rax, %rbx
       add            %rbx, %rcx
       vbroadcastsd   (W, %rax, 8), cos1
       vbroadcastsd   8(W, %rax, 8), sin1
       vbroadcastsd   (W, %rbx, 8), cos2
       vbroadcastsd   8(W, %rbx, 8), sin2
       vbroadcastsd   (W, %rcx, 8), cos3
       vbroadcastsd   8(W, %rcx, 8), sin3
       mov            N2, %rsi
       shl            k2
       sub            k2, %rsi
       shl            k2
       shl            $2, %rsi
       lea            (X, k2, 8), %rax
       lea            (%rax, N, 8), %rbx
       lea            (%rbx, N, 8), %rcx
       lea            (%rcx, N, 8), %rdx
       shr            $2, k2
       vmovapd        (%rax), er0
       vmovapd        (%rbx), er2
       vmovapd        (%rcx), er1
       vmovapd        (%rdx), er3
       vmovapd        (%rax, %rsi, 8), ei0
       vmovapd        (%rbx, %rsi, 8), ei2
       vmovapd        (%rcx, %rsi, 8), ei1
       vmovapd        (%rdx, %rsi, 8), ei3
       vmulpd         sin1, ei1, tr
       vmulpd         cos1, ei1, ti
       vfmsub231pd    cos1, er1, tr
       vfmadd231pd    sin1, er1, ti
       vmovapd        ti, ei1
       vmovapd        tr, er1
       vmulpd         sin2, ei2, tr
       vmulpd         cos2, ei2, ti
       vfmsub231pd    cos2, er2, tr
       vfmadd231pd    sin2, er2, ti
       vmovapd        ti, ei2
       vmovapd        tr, er2
       vmulpd         sin3, ei3, tr
       vmulpd         cos3, ei3, ti
       vfmsub231pd    cos3, er3, tr
       vfmadd231pd    sin3, er3, ti
       vmovapd        ti, ei3
       vmovapd        tr, er3
       vaddpd         er2, er0, tr0
       vaddpd         ei2, ei0, ti0
       vsubpd         er2, er0, tr2
       vsubpd         ei2, ei0, ti2
       vaddpd         er3, er1, tr1
       vaddpd         ei3, ei1, ti1
       vsubpd         er1, er3, tr3
       vsubpd         ei3, ei1, ti3
       vaddpd         tr1, tr0, er0
       vaddpd         ti1, ti0, ei0
       vaddpd         ti3, tr2, er1
       vaddpd         tr3, ti2, ei1
       vsubpd         tr1, tr0, er2
       vsubpd         ti0, ti1, ei2
       vsubpd         ti3, tr2, er3
       vsubpd         ti2, tr3, ei3
       vmovapd        er0, (%rax)
       vmovapd        er1, (%rbx)
       vmovapd        ei3, (%rdx)
       vmovapd        ei2, (%rcx)
       vmovapd        er2, (%rbx, %rsi, 8)
       vmovapd        er3, (%rax, %rsi, 8)
       vmovapd        ei1, (%rcx, %rsi, 8)
       vmovapd        ei0, (%rdx, %rsi, 8)
       jmp            k2loopb
k2loope:
       ret

finish:
       mov            N, NHI
       shr            $2, NHI
       xor            %rcx, %rcx
simd_scramble0_b:
       cmp            NHI, %rcx
       jge            simd_scramble0_e
       mov            SIMD_SIZE, %rax
       mul            %rcx
       mov            %rcx, %rbx
       mov            (X, %rax, 8), %rdx
       mov            %rdx, (Y, %rbx, 8)
       add            $2, %rax
       add            NHI, %rbx
       mov            (X, %rax, 8), %rdx
       mov            %rdx, (Y, %rbx, 8)
       dec            %rax
       add            NHI, %rbx
       mov            (X, %rax, 8), %rdx
       mov            %rdx, (Y, %rbx, 8)
       add            $2, %rax
       add            NHI, %rbx
       mov            (X, %rax, 8), %rdx
       mov            %rdx, (Y, %rbx, 8)
       inc            %rcx
       jmp            simd_scramble0_b
simd_scramble0_e:
       lea            (Y), %rax
       lea            (%rax, N2, 8), %rbx
       lea            (%rbx, N2, 8), %rcx
       lea            (%rcx, N2, 8), %rdx
       vmovq          (%rax), xer0
       vmovq          (%rbx), xer2
       vmovq          (%rcx), xer1
       vmovq          (%rdx), xer3
       vaddpd         xer2, xer0, xtr0
       vsubpd         xer2, xer0, xtr1
       vaddpd         xer3, xer1, xtr2
       vsubpd         xer1, xer3, xtr3
       vaddpd         xtr2, xtr0, xer0
       vsubpd         xtr2, xtr0, xer1
       vmovq          xer0, (%rax)
       vmovq          xtr1, (%rbx)
       vmovq          xer1, (%rcx)
       vmovq          xtr3, (%rdx)
       cmp            $0, N2
       jne            skipr1
       ret

skipr1:
       shr            N2
       lea            (Y, N2, 8), %rax
       shl            N2
       lea            (%rax, N2, 8), %rbx
       lea            (%rbx, N2, 8), %rcx
       lea            (%rcx, N2, 8), %rdx
       vmovq          (%rax), xti0
       vmovq          (%rbx), xti3
       vmovq          (%rcx), xti1
       vmovq          (%rdx), xti2
       vsubpd         xti2, xti1, xer1
       vaddpd         xti2, xti1, xer2
       vmulpd         C, xer1, xti1
       vmulpd         nC, xer2, xti2
       vaddpd         xti1, xti0, xei0
       vsubpd         xti3, xti2, xei1
       vsubpd         xti1, xti0, xei2
       vaddpd         xti2, xti3, xei3
       vmovq          xei1, (%rdx)
       vmovq          xei3, (%rcx)
       vmovq          xei2, (%rbx)
       vmovq          xei0, (%rax)
       xor            k2, k2
       mov            $8, %rdi
       cmp            %rdi, N2
       cmovl          N2, %rdi
       shr            %rdi
k2loopb2:
       inc            k2
       cmp            k2, %rdi
       jle            k2loope2
       mov            k2, %rax
       shl            %rax
       mov            %rax, %rcx
       mov            %rax, %rbx
       add            %rax, %rbx
       add            %rbx, %rcx
       vmovq          (W, %rax, 8), xcos1
       vmovq          8(W, %rax, 8), xsin1
       vmovq          (W, %rbx, 8), xcos2
       vmovq          8(W, %rbx, 8), xsin2
       vmovq          (W, %rcx, 8), xcos3
       vmovq          8(W, %rcx, 8), xsin3
       mov            N2, %rsi
       sub            k2, %rsi
       sub            k2, %rsi
       lea            (Y, k2, 8), %rax
       lea            (%rax, N2, 8), %rbx
       lea            (%rbx, N2, 8), %rcx
       lea            (%rcx, N2, 8), %rdx
       vmovq          (%rax), xer0
       vmovq          (%rbx), xer2
       vmovq          (%rcx), xer1
       vmovq          (%rdx), xer3
       vmovq          (%rax, %rsi, 8), xei0
       vmovq          (%rbx, %rsi, 8), xei2
       vmovq          (%rcx, %rsi, 8), xei1
       vmovq          (%rdx, %rsi, 8), xei3
       vmulpd         xsin1, xei1, xtr
       vmulpd         xcos1, xei1, xti
       vfmsub231pd    xcos1, xer1, xtr
       vfmadd231pd    xsin1, xer1, xti
       vmovq          xti, xei1
       vmovq          xtr, xer1
       vmulpd         xsin2, xei2, xtr
       vmulpd         xcos2, xei2, xti
       vfmsub231pd    xcos2, xer2, xtr
       vfmadd231pd    xsin2, xer2, xti
       vmovq          xti, xei2
       vmovq          xtr, xer2
       vmulpd         xsin3, xei3, xtr
       vmulpd         xcos3, xei3, xti
       vfmsub231pd    xcos3, xer3, xtr
       vfmadd231pd    xsin3, xer3, xti
       vmovq          xti, xei3
       vmovq          xtr, xer3
       vaddpd         xer2, xer0, xtr0
       vaddpd         xei2, xei0, xti0
       vsubpd         xer2, xer0, xtr2
       vsubpd         xei2, xei0, xti2
       vaddpd         xer3, xer1, xtr1
       vaddpd         xei3, xei1, xti1
       vsubpd         xer1, xer3, xtr3
       vsubpd         xei3, xei1, xti3
       vaddpd         xtr1, xtr0, xer0
       vaddpd         xti1, xti0, xei0
       vaddpd         xti3, xtr2, xer1
       vaddpd         xtr3, xti2, xei1
       vsubpd         xtr1, xtr0, xer2
       vsubpd         xti0, xti1, xei2
       vsubpd         xti3, xtr2, xer3
       vsubpd         xti2, xtr3, xei3
       movq           xer0, (%rax)
       movq           xer1, (%rbx)
       movq           xei3, (%rdx)
       movq           xei2, (%rcx)
       movq           xer3, (%rax, %rsi, 8)
       movq           xer2, (%rbx, %rsi, 8)
       movq           xei1, (%rcx, %rsi, 8)
       movq           xei0, (%rdx, %rsi, 8)
       jmp            k2loopb2
k2loope2:
       mov            N, N2
       shr            $2, N2
       xor            k2, k2
       mov            N2, %rdi
       shr            %rdi
k2loopb3:
       add            SIMD_SIZE, k2
       cmp            k2, %rdi
       jle            k2loope3
       mov            k2, %rax
       shl            %rax
       vmovq          %rax, %xmm10
       vpbroadcastq   %xmm10, %ymm10
       vpaddq         iota, %ymm10, %ymm10
       vmovdqa        ones, %ymm13
       vmovdqa        ones, %ymm14
	   vgatherqpd     %ymm13, (W, %ymm10, 8), cos1
	   vgatherqpd     %ymm14, 8(W, %ymm10, 8), sin1
       vmulpd         sin1, sin1, cos2
       vmulpd         cos1, sin1, sin2
       vfmsub231pd    cos1, cos1, cos2
       vfmadd231pd    sin1, cos1, sin2
       vmulpd         sin1, sin2, cos3
       vmulpd         cos1, sin2, sin3
       vfmsub231pd    cos1, cos2, cos3
       vfmadd231pd    sin1, cos2, sin3
       mov            N2, %rsi
       sub            k2, %rsi
       sub            k2, %rsi
       sub            $3, %rsi
       lea            (Y, k2, 8), %rax
       lea            (%rax, N2, 8), %rbx
       lea            (%rbx, N2, 8), %rcx
       lea            (%rcx, N2, 8), %rdx
       vmovapd        (%rax), er0
       vmovapd        (%rbx), er2
       vmovapd        (%rcx), er1
       vmovapd        (%rdx), er3
       vmovdqu        (%rax, %rsi, 8), ei0
       vmovdqu        (%rbx, %rsi, 8), ei2
       vmovdqu        (%rcx, %rsi, 8), ei1
       vmovdqu        (%rdx, %rsi, 8), ei3
       vpermpd        pmask, ei0, ei0
       vpermpd        pmask, ei2, ei2
       vpermpd        pmask, ei1, ei1
       vpermpd        pmask, ei3, ei3
       vmulpd         sin1, ei1, tr
       vmulpd         cos1, ei1, ti
       vfmsub231pd    cos1, er1, tr
       vfmadd231pd    sin1, er1, ti
       vmovapd        ti, ei1
       vmovapd        tr, er1
       vmulpd         sin2, ei2, tr
       vmulpd         cos2, ei2, ti
       vfmsub231pd    cos2, er2, tr
       vfmadd231pd    sin2, er2, ti
       vmovapd        ti, ei2
       vmovapd        tr, er2
       vmulpd         sin3, ei3, tr
       vmulpd         cos3, ei3, ti
       vfmsub231pd    cos3, er3, tr
       vfmadd231pd    sin3, er3, ti
       vmovapd        ti, ei3
       vmovapd        tr, er3
       vaddpd         er2, er0, tr0
       vaddpd         ei2, ei0, ti0
       vsubpd         er2, er0, tr2
       vsubpd         ei2, ei0, ti2
       vaddpd         er3, er1, tr1
       vaddpd         ei3, ei1, ti1
       vsubpd         er1, er3, tr3
       vsubpd         ei3, ei1, ti3
       vaddpd         tr1, tr0, er0
       vaddpd         ti1, ti0, ei0
       vaddpd         ti3, tr2, er1
       vaddpd         tr3, ti2, ei1
       vsubpd         tr1, tr0, er2
       vsubpd         ti0, ti1, ei2
       vsubpd         ti3, tr2, er3
       vsubpd         ti2, tr3, ei3
       vmovapd        er0, (%rax)
       vmovapd        er1, (%rbx)
       vmovapd        ei3, (%rdx)
       vmovapd        ei2, (%rcx)
       vpermpd        pmask, er2, er2
       vpermpd        pmask, er3, er3
       vpermpd        pmask, ei1, ei1
       vpermpd        pmask, ei0, ei0
       vmovdqu        er2, (%rbx, %rsi, 8)
       vmovdqu        er3, (%rax, %rsi, 8)
       vmovdqu        ei1, (%rcx, %rsi, 8)
       vmovdqu        ei0, (%rdx, %rsi, 8)
       jmp            k2loopb3
k2loope3:
       ret

scramble_hi:
       push           W
       bsf            N, %rax
       mov            %rax, %rbx
       shr            %rax
       sub            %rax, %rbx
       mov            $1, NLO
       mov            $1, NHI
       mov            %rax, %rcx
       shl            %rcx, NHI
       mov            %rbx, %rcx
       shl            %rcx, NLO
       jz             skip_hi_scramble
       xor            L, L
       xor            J, J
scramble_hi_loop0:
       xor            K, K
scramble_hi_loop2b:
       mov            L, %rax
       mul            NLO
       add            K, %rax
       shl            $5, %rax
       mov            %rax, %rcx
       mov            J, %rax
       mul            NLO
       add            K, %rax
       shl            $5, %rax
       vmovapd        (Y, %rcx), tmp
       vmovapd        tmp, (X, %rax)
       inc            K
       cmp            K, NLO
       jg             scramble_hi_loop2b
scramble_hi_loop2e:
       mov            NHI, K
       shr            K
       mov            NHI, %rax
       dec            %rax
       cmp            L, %rax
       je             scramble_hi_loop1s
scramble_hi_loop1b:
       cmp            K, J
       jl             scramble_hi_loop1e
       sub            K, J
       shr            K
       jmp            scramble_hi_loop1b
scramble_hi_loop1e:
       add            K, J
scramble_hi_loop1s:
       inc            L
       cmp            L, NHI
       jne            scramble_hi_loop0
skip_hi_scramble:
       pop            W
       ret

transpose:
       push           W
       bsf            N, %rax
       mov            %rax, %rbx
       shr            %rax
       sub            %rax, %rbx
       sub            %rax, %rbx
       mov            $1, NHI
       mov            $1, NMID
       mov            %rax, %rcx
       shl            %rcx, NHI
       mov            %rbx, %rcx
       shl            %rcx, NMID
       xor            J, J
transpose_b0:
       xor            K, K
transpose_b1:
       xor            M, M
transpose_b2:
       mov            M, %rax
       mul            NMID
       add            K, %rax
       mul            NHI
       add            J, %rax
       mov            %rax, %rcx
       mov            J, %rax
       mul            NMID
       add            K, %rax
       mul            NHI
       add            M, %rax
       shl            $5, %rax
       shl            $5, %rcx
	   vmovapd        (X, %rax), %ymm0
	   vmovapd        %ymm0, (Y, %rcx)
       inc            M
       cmp            M, NHI
       jne            transpose_b2
transpose_e2:
       inc            K
       cmp            K, NMID
       jne            transpose_b1
transpose_e1:
       inc            J
       cmp            J, NHI
       jne            transpose_b0
transpose_e0:
       pop            W
       ret
