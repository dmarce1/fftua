#define SIMD_SIZE      $4
#define X              %rdi
#define M              %rsi
#define n0             %r8
#define n1             %r9
#define n              %r10
#define m0             %r11
#define m1             %r12
#define m              %r13
#define A              %r15
#define J0             %ymm0
#define J1             %ymm1
#define J2             %ymm2
#define J3             %ymm3
#define vones          %ymm4
#define B0             %ymm5
#define B1             %ymm6
#define B2             %ymm7
#define B3             %ymm8
#define C0             %ymm9
#define C1             %ymm10
#define C2             %ymm11
#define C3             %ymm12


       .global         fft_scramble

       .data
       .align 32
I0:    .quad  0
       .quad  4
       .quad  8
       .quad  12
I1:    .quad  16
       .quad  20
       .quad  24
       .quad  28
I2:    .quad  32
       .quad  36
       .quad  40
       .quad  44
I3:    .quad  48
       .quad  52
       .quad  56
       .quad  60
ones:  .quad  0xffffffffffffffff
       .quad  0xffffffffffffffff
       .quad  0xffffffffffffffff
       .quad  0xffffffffffffffff
A0:    .double 0.0
       .double 0.0
       .double 0.0
       .double 0.0
A1:    .double 0.0
       .double 0.0
       .double 0.0
       .double 0.0
A2:    .double 0.0
       .double 0.0
       .double 0.0
       .double 0.0
A3:    .double 0.0
       .double 0.0
       .double 0.0
       .double 0.0


       .text
fft_scramble:
       emms
       push            %rbx
       push            %r15
       push            %r14
       push            %r13
       push            %r12
       xor             %rax, %rax
       movq            $4, %rbx
       movq            $8, %rcx
       movq            $12, %r14
       vmovdqa         I0, J0
       vmovdqa         I1, J1
       vmovdqa         I2, J2
       vmovdqa         I3, J3
       xor             %rax, %rax
       not             %rax
       movq            %rax, %xmm4
       vpbroadcastq    %xmm4, vones
       lea             A0, A
       xor             n0, n0
loop_n0:
       mov             n0, m0
       mov             n0, n1
       add             SIMD_SIZE, n1
       cmp             n1, M
       cmovl           M, n1
       mov             n0, n
       cmp             n, n1
       jle             skip_loop_n0
loop_n:
       mov             m0, m1
       add             SIMD_SIZE, m1
       cmp             m1, M
       cmovl           M, m1
       mov             n, m
loop_m:
       inc             m
       cmp             m, m1
       jle             skip_loop_m
       mov             M, %rax
       mul             m
       mov             %rax, %rbx
       mov             M, %rax
       mul             n
       add             n, %rbx
       add             m, %rax
       movq            (X, %rax, 8), %rcx
       movq            (X, %rbx, 8), %r14
       movq            %r14, (X, %rax, 8)
       movq            %rcx, (X, %rbx, 8)
       jmp             loop_m
skip_loop_m:
       inc             n
       cmp             n, n1
       jne             loop_n
       mov             n0, m0
loop_m0:
       add             SIMD_SIZE, m0
       cmp             m0, M
       jle             skip_loop_m0
       mov             n0, %rax
       add             $3, %rax
       mul             M
       mov             %rax, %r14
       mov             n0, %rax
       add             $2, %rax
       mul             M
       mov             %rax, %rcx
       mov             n0, %rax
       inc             %rax
       mul             M
       mov             %rax, %rbx
       mov             n0, %rax
       mul             M
       add             m0, %rax
       add             m0, %rbx
       add             m0, %rcx
       add             m0, %r14
       vmovapd         (X, %rax, 8), C0
       vmovapd         (X, %rbx, 8), C1
       vmovapd         (X, %rcx, 8), C2
       vmovapd         (X, %r14, 8), C3
       vmovapd         C0, A0
       vmovapd         C1, A1
       vmovapd         C2, A2
       vmovapd         C3, A3
       push            %rax
       push            %rbx
       push            %rcx
       push            %r14
       vgatherqpd      vones, (A, J0, 8), B0
       vgatherqpd      vones, (A, J1, 8), B1
       vgatherqpd      vones, (A, J2, 8), B2
       vgatherqpd      vones, (A, J3, 8), B3
       mov             m0, %rax
       add             $3, %rax
       mul             M
       mov             %rax, %r14
       mov             m0, %rax
       add             $2, %rax
       mul             M
       mov             %rax, %rcx
       mov             m0, %rax
       inc             %rax
       mul             M
       mov             %rax, %rbx
       mov             m0, %rax
       mul             M
       add             n0, %rax
       add             n0, %rbx
       add             n0, %rcx
       add             n0, %r14
       vmovapd         (X, %rax, 8), C0
       vmovapd         (X, %rbx, 8), C1
       vmovapd         (X, %rcx, 8), C2
       vmovapd         (X, %r14, 8), C3
       vmovapd         C0, A0
       vmovapd         C1, A1
       vmovapd         C2, A2
       vmovapd         C3, A3
       vmovapd         B0, (X, %rax, 8)
       vmovapd         B1, (X, %rbx, 8)
       vmovapd         B2, (X, %rcx, 8)
       vmovapd         B3, (X, %r14, 8)
       vgatherqpd      vones, (A, J0, 8), B0
       vgatherqpd      vones, (A, J1, 8), B1
       vgatherqpd      vones, (A, J2, 8), B2
       vgatherqpd      vones, (A, J3, 8), B3
       pop             %r14
       pop             %rcx
       pop             %rbx
       pop             %rax
       vmovapd         B0, (X, %rax, 8)
       vmovapd         B1, (X, %rbx, 8)
       vmovapd         B2, (X, %rcx, 8)
       vmovapd         B3, (X, %r14, 8)
       jmp             loop_m0
skip_loop_m0:
       add             SIMD_SIZE, n0
       cmp             n0, M
       jne             loop_n0
skip_loop_n0:
       pop             %r12
       pop             %r13
       pop             %r14
       pop             %r15
       pop             %rbx
       emms
       ret


