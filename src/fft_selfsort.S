#define  SIMD_SIZE      $4
#define  N1             $4

         .global        fft_selfsort

         .data
         .align         32
pM_SQRT1_2:
        .double         +0.70710678118654752
        .double         +0.70710678118654752
        .double         +0.70710678118654752
        .double         +0.70710678118654752
nM_SQRT1_2:
        .double         -0.70710678118654752
        .double         -0.70710678118654752
        .double         -0.70710678118654752
        .double         -0.70710678118654752


#define  X              %r8
#define  W              %r9
#define  N              %r10

#define  log2MMM        %r11
#define  log2MMID       %r12
#define  log2MA         %r13

fft_selfsort:
         push           %r12
         push           %r13


loop0:
         push           X
         pop            X
         jnz            loop0

         pop            %r13
         pop            %r12
         ret



fft_swap_indices:
         xor            %rax, %rax
         mov            log2MMM, %rcx
         mov            $1, %rbx
         xor            %rdx, %rdi
         xor            %rbx, %rsi
         shl            %cl, %rbx
swap4:
         push           %rdi
         push           %rsi
         push           %rax
         push           %rbx
         add            %rax, %rdi
         add            %rax, %rsi
         xor            %rax, %rax
         mov            log2MA, %rcx
         mov            $1, %rbx
         shl            %cl, %rdi
         shl            %cl, %rsi
         shl            %cl, %rbx
swap3:
         push           %rdi
         push           %rsi
         push           %rax
         push           %rbx
         add            %rax, %rdi
         add            %rax, %rsi
         xor            %rax, %rax
         mov            log2MMID, %rcx
         mov            $1, %rbx
         shl            %cl, %rdi
         shl            %cl, %rsi
         shl            %cl, %rbx
swap2:
         push           %rdi
         push           %rsi
         push           %rax
         push           %rbx
         add            %rax, %rdi
         add            %rax, %rsi
         xor            %rax, %rax
         mov            log2MA, %rcx
         mov            $1, %rbx
         shl            %cl, %rdi
         shl            %cl, %rsi
         shl            %cl, %rbx
swap1:
         push           %rdi
         push           %rsi
         push           %rax
         push           %rbx
         add            %rax, %rdi
         add            %rax, %rsi
         xor            %rax, %rax
         mov            log2MMM, %rcx
         mov            $1, %rbx
         shl            %cl, %rdi
         shl            %cl, %rsi
         shl            %cl, %rbx
swap0:
         push           %rdi
         push           %rsi
         push           %rax
         push           %rbx
         add            %rax, %rdi
         add            %rax, %rsi
         mov            (X, %rsi, 8), %rbx
         mov            (X, %rdi, 8), %rdx
         mov            %rdx, (X, %rsi, 8)
         mov            %rbx, (X, %rdi, 8)
         pop            %rbx
         pop            %rax
         pop            %rsi
         pop            %rdi
         inc            %rax
         cmp            %rax, %rbx
         jne            swap0
         pop            %rbx
         pop            %rax
         pop            %rsi
         pop            %rdi
         inc            %rax
         cmp            %rax, %rbx
         jne            swap1
         pop            %rbx
         pop            %rax
         pop            %rsi
         pop            %rdi
         inc            %rax
         cmp            %rax, %rbx
         jne            swap2
         pop            %rbx
         pop            %rax
         pop            %rsi
         pop            %rdi
         inc            %rax
         cmp            %rax, %rbx
         jne            swap3
         pop            %rbx
         pop            %rax
         pop            %rsi
         pop            %rdi
         inc            %rax
         cmp            %rax, %rbx
         jne            swap4
         ret

#define  d8N2NLO        %r8
#define  X              %r9
#define  W              %r10
#define  N              %r11
#define  NLO            %r12
#define  nlo            %r13
#define  N2             %r14
#define  k2             %r15
#define  er0            %ymm8
#define  ei0            %ymm9
#define  er1            %ymm10
#define  ei1            %ymm11
#define  er2            %ymm12
#define  ei2            %ymm13
#define  er3            %ymm14
#define  ei3            %ymm15
#define  tr0            %ymm0
#define  ti0            %ymm1
#define  tr1            %ymm2
#define  ti1            %ymm3
#define  tr2            %ymm4
#define  ti2            %ymm5
#define  tr3            %ymm6
#define  ti3            %ymm7
#define  cos1           -8(%ebp)
#define  sin1           -16(%ebp)
#define  cos2           -24(%ebp)
#define  sin2           -32(%ebp)
#define  cos3           -40(%ebp)
#define  sin3           -48(%ebp)
#define  DI             -56(%ebp)
#define  STACK_SIZE     56

fft_apply_butterfly4:
         push           %rbp
         mov            %rsp, %rbx
         add            $STACK_SIZE, %rcx
         mov            %rsp, %rbp
         sub            %rsp, %rcx
         push           %rbx
         push           %r12
         push           %r13
         push           %r14
         push           %r15
         mov            %rdi, X
         mov            %rsi, W
         mov            %rdx, N
         mov            %rcx, NLO
         mov            %r8,  N2
         mov            NLO, %rax
         mul            N2
         shl            $3, %rax
         mov            %rax, d8N2NLO
         xor            nlo, nlo
nloloop0:
         mov            nlo, %rax
         lea            (X, %rax, 8), %rax
         lea            (%rax, d8N2NLO), %rbx
         lea            (%rbx, d8N2NLO), %rcx
         lea            (%rcx, d8N2NLO), %rdx
         vmovapd        (%rax), er0
         vmovapd        (%rbx), er1
         vmovapd        (%rcx), er2
         vmovapd        (%rdx), er3
         vaddpd         er1, er0, tr0
         vsubpd         er1, er0, tr1
         vaddpd         er3, er2, tr2
         vsubpd         er3, er2, tr3
         vmovapd        er0, (%rax)
         vmovapd        tr1, (%rbx)
         vmovapd        er1, (%rcx)
         vmovapd        tr3, (%rdx)
         add            SIMD_SIZE, nlo
         cmp            NLO, nlo
         jne            nloloop0
         cmp            $1, N2
         jle            skipN2o2
         xor            nlo, nlo
nloloop1:
         mov            N2, %rax
         shr            %rax
         add            nlo, %rax
         lea            (X, %rax, 8), %rax
         lea            (%rax, d8N2NLO), %rbx
         lea            (%rbx, d8N2NLO), %rcx
         lea            (%rcx, d8N2NLO), %rdx
         vmovapd        (%rax), tr0
         vmovapd        (%rbx), tr3
         vmovapd        (%rcx), tr1
         vmovapd        (%rdx), tr2
         vsubpd         tr2, tr1, er1
         vaddpd         tr2, tr1, er2
         vmulpd         pM_SQRT1_2, er1, tr1
         vmulpd         nM_SQRT1_2, er2, tr2
         vaddpd         tr1, tr0, er0
         vsubpd         tr3, tr2, er1
         vsubpd         tr1, tr0, er2
         vaddpd         tr2, tr3, er3
         vmovapd        er1, (%rdx)
         vmovapd        er3, (%rcx)
         vmovapd        er2, (%rbx)
         vmovapd        er0, (%rax)
         add            SIMD_SIZE, nlo
         cmp            NLO, nlo
         jne            nloloop1
skipN2o2:
         mov            $1, k2
k2loop:
         mov            N2, %rax
         sub            k2, %rax
         sub            k2, %rax
         mul            NLO
         shl            $3, %rax
         mov            %rax, DI
         mov            k2, %rax
         shl            %rax
         mov            %rax, %rbx
         mov            %rax, %rcx
         add            %rax, %rbx
         add            %rbx, %rcx
         mov            (W, %rax, 8), %rdx
         mov            8(W, %rax, 8), %rax
         mov            %rdx, cos1
         mov            %rax, sin1
         mov            (W, %rbx, 8), %rdx
         mov            8(W, %rbx, 8), %rax
         mov            %rdx, cos2
         mov            %rax, sin2
         mov            (W, %rcx, 8), %rdx
         mov            8(W, %rcx, 8), %rax
         mov            %rdx, cos3
         mov            %rax, sin3
         xor            nlo, nlo
nloloop2:
         vbroadcastsd   cos1, tr1
         vbroadcastsd   sin1, ti1
         vbroadcastsd   cos2, tr2
         vbroadcastsd   sin2, ti2
         vbroadcastsd   cos3, tr3
         vbroadcastsd   sin3, ti3
         mov            k2, %rax
         mul            NLO
         add            nlo, %rax
         lea            (X, %rax, 8), %rax
         lea            (%rax, d8N2NLO), %rbx
         lea            (%rbx, d8N2NLO), %rcx
         lea            (%rcx, d8N2NLO), %rdx
         vmovapd        (%rax), er0
         vmovapd        (%rbx), er2
         vmovapd        (%rcx), er1
         vmovapd        (%rdx), er3
         mov            DI, %rdi
         vmovapd        (%rax, %rdi), ei0
         vmovapd        (%rbx, %rdi), ei2
         vmovapd        (%rcx, %rdi), ei1
         vmovapd        (%rdx, %rdi), ei3
         vmulpd         tr1, ei1, tr0
         vmulpd         ti1, ei1, ti0
         vfmsub231pd    ti1, er1, tr0
         vfmadd231pd    tr1, er1, ti0
         vmovapd        ti0, ei1
         vmovapd        tr0, er1
         vmulpd         tr2, ei2, tr0
         vmulpd         ti2, ei2, ti0
         vfmsub231pd    ti2, er2, tr0
         vfmadd231pd    tr2, er2, ti0
         vmovapd        ti0, ei2
         vmovapd        tr0, er2
         vmulpd         tr3, ei3, tr0
         vmulpd         ti3, ei3, ti0
         vfmsub231pd    ti3, er3, tr0
         vfmadd231pd    tr3, er3, ti0
         vmovapd        ti0, ei3
         vmovapd        tr0, er3
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
         vmovapd        er2, (%rbx, %rdi)
         vmovapd        er3, (%rax, %rdi)
         vmovapd        ei0, (%rdx, %rdi)
         vmovapd        ei1, (%rcx, %rdi)
         vmovapd        ei2, (%rcx)
         vmovapd        ei3, (%rdx)
         add            SIMD_SIZE, nlo
         cmp            NLO, nlo
         jne            nloloop2
         inc            k2
         mov            N2, %rax
         shr            %rax
         cmp            k2, %rax
         jne            k2loop
         pop            %r15
         pop            %r14
         pop            %r13
         pop            %r12
         pop            %rbx
         pop            %rax
         mov            %rbp, %rsp
         pop            %rbp
         ret

mask_inc:
         mov            %rsi, %rax
         not            %rax
         or             %rdi, %rax
         inc            %rax
         and            %rsi, %rax
         ret

