#define  STACK_SIZE  $32
#define  In          %r8
#define  Ia          %r9
#define  Ib          %r10
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
#define  log2NA2     -25(%rbp)
#define  shiftA      -26(%rbp)
#define  shiftB      -27(%rbp)

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
         mov       NA, %rax
         mul       NA
         bsf       %rax, %rax
         mov       %al, log2NA2
         mov       NLO, %rax
         mov       %rax, %rbx
         dec       %rbx
         not       %rbx
         mul       NA
         mov       %rax, %rcx
         dec       %rcx
         and       %rbx, %rcx
         mov       %rcx, %rdx
         mov       NMID, %rdi
         mul       %rdi
         mov       %rax, %rbx
         dec       %rbx
         not       %rbx
         mul       NA
         mov       %rax, %rcx
         dec       %rcx
         and       %rbx, %rcx
         or        %rcx, %rdx
         mov       %rdx, mask
         not       %rdx
         mov       %rdx, not_mask
         bsf       NLO, %rdx
         mov       %rdx, %rcx
         bsf       NA, %rax
         add       %rax, %rdx
         bsf       NMID, %rax
         add       %rax, %rdx
         mov       %cl, shiftA
         mov       %dl, shiftB
         xor       In, In
loop_hi:
         xor       Ia, Ia
loop_a:
         mov       Ia, %rax
         mov       shiftA, %cl
         shl       %cl, %rax
         mov       Ia, %rbx
         mov       shiftB, %cl
         shl       %cl, %rbx
         xor       Ib, Ib
loop_bb:
         cmp       Ib, Ia
         jle       loop_be
         mov       Ib, %rdx
         mov       shiftB, %cl
         shl       %cl, %rdx
         or        %rdx, %rax
         mov       Ib, %rdx
         mov       shiftA, %cl
         shl       %cl, %rdx
         or        %rdx, %rbx
         or        In, %rax
         or        In, %rbx
         mov       (X, %rax, 8), %rcx
         mov       (X, %rbx, 8), %rdx
         mov       %rdx, (X, %rax, 8)
         mov       %rcx, (X, %rbx, 8)
         inc       Ib
         jmp       loop_bb
loop_be:
         inc       Ia
         cmp       Ia, NA
         jne       loop_a
         mov       In, %rax
         or        mask, %rax
         inc       %rax
         and       not_mask, %rax
         mov       %rax, In
         cmp       %rax, N
         jg        loop_hi
         pop       %r15
         pop       %r14
         pop       %r13
         pop       %r12
         pop       %rbx
         leave
         ret
