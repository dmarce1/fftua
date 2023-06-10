#define  SIMD_SIZE   $4
#define  STACK_SIZE  $32
#define  Ia          %r8
#define  Ib          %r9
#define  In          %r10
#define  X           %r11
#define  N           %r12
#define  NA          %r13
#define  mask        %r14
#define  not_mask    %r15
#define  Ja          %rdi
#define  Jb          %rsi
#define  NHI         -8(%rbp)
#define  NMID        -16(%rbp)
#define  NLO         -24(%rbp)
#define  shiftA      -25(%rbp)
#define  shiftB      -26(%rbp)

         .global   fft_swap

fft_swap:
         enter     STACK_SIZE, $0
         push      %rbx
         push      %r12
         push      %r13
         push      %r14
         push      %r15
         mov       %rdi, X
         mov       %rsi, N
         mov       %rdx, NHI
         mov       %rcx, NMID
         mov       %r8, NLO
         mov       %r9, NA
         mov       NLO, %rax
         mov       %rax, %rbx
         dec       %rbx
         not       %rbx
         mul       NA
         mov       %rax, %rcx
         dec       %rcx
         and       %rbx, %rcx
         mov       %rcx, %rsi
         mov       NMID, %rdi
         mul       %rdi
         mov       %rax, %rbx
         dec       %rbx
         not       %rbx
         mul       NA
         mov       %rax, %rcx
         dec       %rcx
         and       %rbx, %rcx
         or        %rcx, %rsi
         mov       %rsi, mask
         not       %rsi
         mov       %rsi, not_mask
         bsf       NLO, %rdx
         mov       %rdx, %rcx
         bsf       NA, %rax
         add       %rax, %rdx
         bsf       NMID, %rax
         add       %rax, %rdx
         mov       %cl, shiftA
         mov       %dl, shiftB
         xor       In, In
loop_hi2:
         xor       Ia, Ia
loop_a2:
         xor       Ib, Ib
loop_bb2:
         cmp       Ib, Ia
         jle       loop_be2
         mov       Ia, %rax
         mov       Ib, %rdx
         mov       shiftA, %cl
         shl       %cl, %rax
         shl       %cl, %rdx
         or        %rdx, %rax
         or        In, %rax
         mov       Ib, %rbx
         mov       Ia, %rdx
         mov       shiftB, %cl
         shl       %cl, %rbx
         shl       %cl, %rdx
         or        %rdx, %rbx
         or        In, %rbx
         mov       (X, %rax, 8), %rcx
         mov       (X, %rbx, 8), %rdx
         mov       %rdx, (X, %rax, 8)
         mov       %rcx, (X, %rbx, 8)
         inc       Ib
         jmp       loop_bb2
loop_be2:
         inc       Ia
         cmp       Ia, NA
         jg        loop_a2
         mov       In, %rax
         or        mask, %rax
         inc       %rax
         and       not_mask, %rax
         mov       %rax, In
         cmp       %rax, N
         jg        loop_hi2
         pop       %r15
         pop       %r14
         pop       %r13
         pop       %r12
         pop       %rbx
         leave
         ret
