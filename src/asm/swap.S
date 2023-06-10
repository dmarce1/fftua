#define  SIMD_SIZE   $4
#define  STACK_SIZE  $72
#define  X           %r8
#define  In          %r9
#define  Ia1         %r10
#define  Ib1         %r11
#define  Ia2         %r12
#define  Ib2         %r13
#define  mask        %r14
#define  not_mask    %r15
#define  NHI         -8(%rbp)
#define  NMID        -16(%rbp)
#define  NLO         -24(%rbp)
#define  NA          -32(%rbp)
#define  N           -40(%rbp)
#define  maxA        -48(%rbp)
#define  maxB        -56(%rbp)
#define  dA          -64(%rbp)
#define  dB          -72(%rbp)

         .global   fft_swap

fft_swap:
         enter     STACK_SIZE, $0
         push      %rbx
         push      %r12
         push      %r13
         push      %r14
         push      %r15
         mov       %r8, NLO
         mov       %r9, NA
         mov       %rdi, X
         mov       %rsi, N
         mov       %rdx, NHI
         mov       %rcx, NMID
         mov       NLO, %rax
         mov       %rax, dA
         mov       %rax, %rbx
         dec       %rbx
         not       %rbx
         mov       NA, %rdi
         mul       %rdi
         mov       %rax, maxA
         mov       %rax, %rcx
         dec       %rcx
         and       %rbx, %rcx
         mov       %rcx, %rsi
         mov       NMID, %rdi
         mul       %rdi
         mov       %rax, dB
         mov       %rax, %rbx
         dec       %rbx
         not       %rbx
         mov       NA, %rdi
         mul       %rdi
         mov       %rax, maxB
         mov       %rax, %rcx
         dec       %rcx
         and       %rbx, %rcx
         or        %rcx, %rsi
         mov       %rsi, mask
         not       %rsi
         mov       %rsi, not_mask
         xor       In, In
         cmp       $4, NLO
         jge       simd_b
loop_hi1:
         xor       Ia1, Ia1
         xor       Ia2, Ia2
loop_a1:
         xor       Ib1, Ib1
         xor       Ib2, Ib2
loop_bb1:
         cmp       Ib1, Ia1
         jle       loop_be2
         mov       In, %rax
         mov       In, %rbx
         or        Ia1, %rax
         or        Ib1, %rax
         or        Ia2, %rbx
         or        Ib2, %rbx
         mov       (X, %rax, 8), %rcx
         mov       (X, %rbx, 8), %rdx
         mov       %rdx, (X, %rax, 8)
         mov       %rcx, (X, %rbx, 8)
         add       dB, Ib1
         add       dA, Ib2
         jmp       loop_bb1
loop_be1:
         add       dA, Ia1
         add       dB, Ia2
         cmp       Ia1, maxA
         jg        loop_a1
         mov       In, %rax
         or        mask, %rax
         inc       %rax
         and       not_mask, %rax
         mov       %rax, In
         cmp       %rax, N
         jg        loop_hi1
         jmp       simd_e
simd_b:
loop_hi2:
         xor       Ia1, Ia1
         xor       Ia2, Ia2
loop_a2:
         xor       Ib1, Ib1
         xor       Ib2, Ib2
loop_bb2:
         cmp       Ib1, Ia1
         jle       loop_be2
         mov       In, %rax
         mov       In, %rbx
         or        Ia1, %rax
         or        Ib1, %rax
         or        Ia2, %rbx
         or        Ib2, %rbx
         vmovapd   (X, %rax, 8), %ymm0
         vmovapd   (X, %rbx, 8), %ymm1
         vmovapd   %ymm1, (X, %rax, 8)
         vmovapd   %ymm0, (X, %rbx, 8)
         add       dB, Ib1
         add       dA, Ib2
         jmp       loop_bb2
loop_be2:
         add       dA, Ia1
         add       dB, Ia2
         cmp       Ia1, maxA
         jg        loop_a2
         mov       In, %rax
         or        mask, %rax
         add       SIMD_SIZE, %rax
         and       not_mask, %rax
         mov       %rax, In
         cmp       %rax, N
         jg        loop_hi2
simd_e:
         pop       %r15
         pop       %r14
         pop       %r13
         pop       %r12
         pop       %rbx
         leave
         ret
