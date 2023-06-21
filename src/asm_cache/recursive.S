#define    M_SQRT1_2     0.70710678118654752440
#define    N1            $4
#define    Ntw           %r8
#define    m             %r9
#define    M             %r10
#define    k2            %r11
#define    N             %r12
#define    X             %r13
#define    C             %r14
#define    S             %r15
#define    er0           %ymm0
#define    ei0           %ymm1
#define    er1           %ymm2
#define    ei1           %ymm3
#define    er2           %ymm4
#define    ei2           %ymm5
#define    er3           %ymm6
#define    ei3           %ymm7
#define    tr0           %ymm8
#define    tr1           %ymm9
#define    tr2           %ymm10
#define    tr3           %ymm11
#define    ti0           %ymm12
#define    ti1           %ymm13
#define    ti2           %ymm14
#define    ti3           %ymm15
#define    cos1          -32(%rbp)
#define    cos2          -64(%rbp)
#define    cos3          -128(%rbp)
#define    sin1          -160(%rbp)
#define    sin2          -192(%rbp)
#define    sin3          -224(%rbp)
#define    STACK_SIZE    $224


           .global       fft_recursive_hi


           .data

           .align        32
tw45:      .double       M_SQRT1_2
           .double       M_SQRT1_2
           .double       M_SQRT1_2
           .double       M_SQRT1_2
none:      .double       -1.0
           .double       -1.0
           .double       -1.0
           .double       -1.0


           .text

fft_recursive_hi:
           enter         STACK_SIZE, $0
           push          %r15
           push          %r14
           push          %r13
           push          %r12
           push          %rbx
           mov           %r8, M
           mov           %rcx, N
           mov           %rdx, S
           mov           %rsi, C
           mov           %rdi, X
           mov           $1, Ntw
           call          recurse
           pop           %rbx
           pop           %r12
           pop           %r13
           pop           %r14
           pop           %r15
           leave
           ret

recurse:   cmp           $2, N
           jne           rcalls
           xor           m, m
r2_loop:   lea           (X, m, 8), %rax
           mov           M, %rdx
           imul          N, %rdx
           lea           (%rax, %rdx, 4), %rbx
           vmovapd       (%rax), er0
           vmovapd       (%rbx), er1
           vaddpd        er1, er0, tr0
           vsubpd        er1, er0, tr1
           vmovapd       tr0, (%rax)
           vmovapd       tr1, (%rbx)
           add           $4, m
           cmp           m, M
           jg            r2_loop
           jmp           rdone
rcalls:    cmp           $4, N
           je            skip
           mov           N, %rax
           imul          M, %rax
           push          X
           lea           (X, %rax, 2), X
           push          X
           lea           (X, %rax, 2), X
           push          X
           lea           (X, %rax, 2), X
           shr           $2, N
           shl           $2, Ntw
           call          recurse
           pop           X
           call          recurse
           pop           X
           call          recurse
           pop           X
           call          recurse
           shr           $2, Ntw
           shl           $2, N
skip:      xor           m, m
k0_loop:   lea           (X, m, 8), %rax
           mov           M, %rdx
           imul          N, %rdx
           lea           (%rax, %rdx, 2), %rbx
           lea           (%rbx, %rdx, 2), %rcx
           lea           (%rcx, %rdx, 2), %rdx
           vmovapd       (%rax), er0
           vmovapd       (%rbx), er2
           vmovapd       (%rcx), er1
           vmovapd       (%rdx), er3
           vaddpd        er2, er0, tr0
           vsubpd        er2, er0, tr2
           vaddpd        er3, er1, tr1
           vsubpd        er1, er3, tr3
           vaddpd        tr1, tr0, er0
           vsubpd        tr1, tr0, er2
           vmovapd       er0, (%rax)
           vmovapd       tr2, (%rbx)
           vmovapd       er2, (%rcx)
           vmovapd       tr3, (%rdx)
           add           $4, m
           cmp           m, M
           jg            k0_loop
           cmp           $8, N
           jl            rdone
           xor           m, m
Ny_loop:   mov           M, %rdx
           imul          N, %rdx
           lea           (X, %rdx), %rax
           lea           (%rax, m, 8), %rax
           lea           (%rax, %rdx, 2), %rbx
           lea           (%rbx, %rdx, 2), %rcx
           lea           (%rcx, %rdx, 2), %rdx
           vmovapd       (%rax), er0
           vmovapd       (%rbx), er2
           vmovapd       (%rcx), er1
           vmovapd       (%rdx), er3
           vaddpd        er3, er1, tr0
           vsubpd        er3, er1, tr2
           vmulpd        tw45, tr2, tr1
           vmulpd        tw45, tr0, tr3
           vmovapd       er0, tr0
           vmovapd       er2, tr2
           vaddpd        tr1, tr0, er0
           vaddpd        tr3, tr2, er3
           vsubpd        tr1, tr0, er1
           vsubpd        tr3, tr2, er2
           vmulpd        none, er3, er3
           vmovapd       er0, (%rax)
           vmovapd       er1, (%rbx)
           vmovapd       er2, (%rcx)
           vmovapd       er3, (%rdx)
           add           $4, m
           cmp           m, M
           jg            Ny_loop
           mov           M, %rdi
           xor           k2, k2
k2_loop_b: inc           k2
           mov           N, %rax
           shr           $3, %rax
           cmp           k2, %rax
           jle           rdone
           mov           N, %rsi
           shr           $2, %rsi
           sub           k2, %rsi
           sub           k2, %rsi
           imul          M, %rsi
           mov           k2, %rax
           imul          Ntw, %rax
           imul          $2, %rax, %rbx
           imul          $3, %rax, %rcx
           push          %rdi
           push          %rsi
           movq          (C, %rax, 8), %rdx
           movq          (C, %rbx, 8), %rdi
           movq          (C, %rcx, 8), %rsi
           movq          (S, %rax, 8), %rax
           movq          (S, %rbx, 8), %rbx
           movq          (S, %rcx, 8), %rcx
           movq          %rdx, cos1
           movq          %rdi, cos2
           movq          %rsi, cos3
           movq          %rax, sin1
           movq          %rbx, sin2
           movq          %rcx, sin3
           pop           %rsi
           pop           %rdi
           xor           m, m
m_loop:    vbroadcastsd  cos1, tr1
           vbroadcastsd  cos2, tr2
           vbroadcastsd  cos3, tr3
           vbroadcastsd  sin1, ti1
           vbroadcastsd  sin2, ti2
           vbroadcastsd  sin3, ti3
           mov           M, %rdx
           imul          N, %rdx
           lea           (X, %rdi, 8), %rax
           lea           (%rax, %rdx, 2), %rbx
           lea           (%rbx, %rdx, 2), %rcx
           lea           (%rcx, %rdx, 2), %rdx
           vmovapd       (%rax), er0
           vmovapd       (%rbx), er2
           vmovapd       (%rcx), er1
           vmovapd       (%rdx), er3
           vmovapd       (%rax, %rsi, 8), ei0
           vmovapd       (%rbx, %rsi, 8), ei2
           vmovapd       (%rcx, %rsi, 8), ei1
           vmovapd       (%rdx, %rsi, 8), ei3
           vmulpd        ti1, ei1, tr0
           vmulpd        tr1, ei1, ti0
           vfmsub231pd   tr1, er1, tr0
           vfmadd231pd   ti1, er1, ti0
           vmovapd       ti0, ei1
           vmovapd       tr0, er1
           vmulpd        ti2, ei2, tr0
           vmulpd        tr2, ei2, ti0
           vfmsub231pd   tr2, er2, tr0
           vfmadd231pd   ti2, er2, ti0
           vmovapd       ti0, ei2
           vmovapd       tr0, er2
           vmulpd        ti3, ei3, tr0
           vmulpd        tr3, ei3, ti0
           vfmsub231pd   tr3, er3, tr0
           vfmadd231pd   ti3, er3, ti0
           vmovapd       ti0, ei3
           vmovapd       tr0, er3
           vaddpd        er2, er0, tr0
           vaddpd        ei2, ei0, ti0
           vsubpd        er2, er0, tr2
           vsubpd        ei2, ei0, ti2
           vaddpd        er3, er1, tr1
           vaddpd        ei3, ei1, ti1
           vsubpd        er1, er3, tr3
           vsubpd        ei3, ei1, ti3
           vaddpd        tr1, tr0, er0
           vaddpd        ti1, ti0, ei0
           vaddpd        ti3, tr2, er1
           vaddpd        tr3, ti2, ei1
           vsubpd        tr1, tr0, er2
           vsubpd        ti0, ti1, ei2
           vsubpd        ti3, tr2, er3
           vsubpd        ti2, tr3, ei3
           vmovapd       er0, (%rax)
           vmovapd       er1, (%rbx)
           vmovapd       ei2, (%rcx)
           vmovapd       ei3, (%rdx)
           vmovapd       er3, (%rax, %rsi, 8)
           vmovapd       er2, (%rbx, %rsi, 8)
           vmovapd       ei1, (%rcx, %rsi, 8)
           vmovapd       ei0, (%rdx, %rsi, 8)
           add           $4, %rdi
           add           $4, m
           cmp           m, M
           jg            m_loop
           jmp           k2_loop_b
rdone:     ret

