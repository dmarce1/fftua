#define  STACK_SIZE   $8
#define  N1           $4
#define  X            %r12
#define  W            %r13
#define  N2           %r15
#define  NMID         %r14
#define  NHI          %rbx
#define  N            -8(%rbp)

         .global      fft_selfsort

         .text

fft_selfsort:
         enter        STACK_SIZE, $0
         push         %rbx
         push         %r12
         push         %r13
         push         %r14
         push         %r15
         mov          %rdi, X
         mov          %rsi, W
         mov          %rdx, N
         mov          $1, N2
         mov          N2, NHI
         mov          X, %rdi
         mov          N, %rsi
         call         trivial_butterfly
         shl          $2, NHI
         shl          $2, N2
         mov          N2, %rdx
         mov          NHI, %rax
         mul          %rdx
         mov          %rax, %rcx
         shl          $4, %rcx
         mov          N, %rax
         xor          %rdx, %rdx
         div          %rcx
         mov          %rax, NMID
L0:
         cmp          $0, NMID
         je           L2
         mov          X, %rdi
         mov          W, %rsi
         mov          NHI, %rdx
         mov          NMID, %rcx
         mov          N2, %r8
         call         butterfly4_and_transpose
         shl          $2, NHI
         shl          $2, N2
         mov          N2, %rdx
         mov          NHI, %rax
         mul          %rdx
         mov          %rax, %rcx
         shl          $4, %rcx
         mov          N, %rax
         xor          %rdx, %rdx
         div          %rcx
         mov          %rax, NMID
         jmp          L0
L2:
         mov          N2, %rax
         mov          N, %rcx
         shl          $2, %rax
         shr          $2, %rcx
         cmp          %rax, %rcx
         jl           L3
         mov          N, %rax
         xor          %rdx, %rdx
         mov          N2, %rcx
         shl          $2, %rcx
         div          %rcx
         mov          %rax, NHI
         mov          X, %rdi
         mov          W, %rsi
         mov          NHI, %rdx
         mov          N2, %rcx
	     call         butterfly4
         shl          $2, N2
         jmp          L2
L3:
         mov          X, %rdi
         mov          N, %rsi
         shr          $4, %rsi
         call         transpose
         mov          N, %rax
         xor          %rdx, %rdx
         mov          N2, %rcx
         shl          $2, %rcx
         div          %rcx
         mov          %rax, NHI
         mov          X, %rdi
         mov          W, %rsi
         mov          N, %rdx
         mov          N2, %rcx
         call         butterfly4_finish
         pop          %r15
         pop          %r14
         pop          %r13
         pop          %r12
         pop          %rbx
         leave
         ret
