#define M_SQRT2_1     0.70710678118
#define N1               $4
#define SIMD_SIZE     $4
#define X             %r8
#define Y             %r9
#define W             %r10
//#define Wv            %r11
#define NHI           %r12
#define N             %r13
#define N2            %r14
#define k2            %r15
#define I             %rax
#define J             %rbx
#define I4            %rcx
#define J4            %rdx
#define K             %rdi
#define Nm1           %rsi
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

       .global        fft_recursive

       .data
       .align         32
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
Z:
Z0:    .double        0.0
Z1:    .double        0.0
Z2:    .double        0.0
Z3:    .double        0.0


       .text

fft_recursive:
       push           %rbx
       push           %r12
       push           %r13
       push           %r14
       push           %r15
       mov            %rdi, Y
       mov            %rsi, X
       mov            %rdx, W
       mov            %rcx, N
       mov            N, N2
       mov            $1, NHI
       shr            $2, N2
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
       mov            N, Nm1
       dec            Nm1
       xor            I, I
       xor            J, J
       xor            I4, I4
       xor            J4, J4
scramble_loop0:
       vmovapd        (Y, I4, 8), tmp
       vmovapd        tmp, (X, J4, 8)
       mov            N, K
       shr            K
scramble_loop1b:
       cmp            K, J
       jl             scramble_loop1e
       sub            K, J
       shr            K
       jmp            scramble_loop1b
scramble_loop1e:
       add            K, J
       mov            J, J4
       shl            $2, J4
       inc            I
       add            SIMD_SIZE, I4
       cmp            I, Nm1
       jne            scramble_loop0
       shl            $2, Nm1
       vmovapd        (Y, Nm1, 8), tmp
       vmovapd        tmp, (X, Nm1, 8)
       call           recurse
       pop            %r15
       pop            %r14
       pop            %r13
       pop            %r12
       pop            %rbx
       emms
       ret
recurse:
       cmp            $4, N
       jne            Nneq4
       jmp            skipcalls
Nneq4:
       cmp            $8, N
       jne            Nneq8
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
       jmp            skipcalls
Nneq8:
       mov            N, %rbx
       shr            $2, N
       shr            $2, N2
       shl            $2, NHI
       push           X
       lea            (X, %rbx, 8), X
       push           X
       lea            (X, %rbx, 8), X
       push           X
       lea            (X, %rbx, 8), X
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
skipcalls:
       shl            $2, N2
       lea            (X), %rax
       lea            (%rax, N2, 8), %rbx
       lea            (%rbx, N2, 8), %rcx
       lea            (%rcx, N2, 8), %rdx
       shr            $2, N2
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
       jle            not_k2eqN2o2
       shl            N2
       lea            (X, N2, 8), %rax
       shl            N2
       lea            (%rax, N2, 8), %rbx
       lea            (%rbx, N2, 8), %rcx
       lea            (%rcx, N2, 8), %rdx
       shr            $2, N2
       vmovapd        (%rax), ti0
       vmovapd        (%rbx), ti3
       vmovapd        (%rcx), ti1
       vmovapd        (%rdx), ti2
       vsubpd         ti2, ti1, er1
       vaddpd         ti2, ti1, er2
       vmulpd         C, er1, ti1
       vmulpd         nC, er2, ti2
       vaddpd         ti1, ti0, ei0
       vsubpd         ti3, ti2, ei1
       vsubpd         ti1, ti0, ei2
       vaddpd         ti2, ti3, ei3
       vmovapd        ei1, (%rdx)
       vmovapd        ei3, (%rcx)
       vmovapd        ei2, (%rbx)
       vmovapd        ei0, (%rax)
not_k2eqN2o2:
       xor            k2, k2
k2loopb:
       mov            N2, %rcx
       inc            k2
       shr            %rcx
       cmp            k2, %rcx
       jle            k2loope
       mov            NHI, %rax
       mul            k2
       shl            %rax
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
       sub            k2, %rsi
       sub            k2, %rsi
       shl            $2, %rsi
       shl            $2, k2
       shl            $2, N2
       lea            (X, k2, 8), %rax
       lea            (%rax, N2, 8), %rbx
       lea            (%rbx, N2, 8), %rcx
       lea            (%rcx, N2, 8), %rdx
       shr            $2, k2
       shr            $2, N2
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
       vmovapd        er2, (%rbx, %rsi, 8)
       vmovapd        er3, (%rax, %rsi, 8)
       vmovapd        ei3, (%rdx)
       vmovapd        ei2, (%rcx)
       vmovapd        ei1, (%rcx, %rsi, 8)
       vmovapd        ei0, (%rdx, %rsi, 8)
       jmp            k2loopb
k2loope:
       ret

