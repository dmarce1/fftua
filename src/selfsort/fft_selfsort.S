#define  STACK_SIZE   $24
#define  N1           $4
#define  X            %r12
#define  C            %r13
#define  S            %r14
#define  N2           %r15
#define  NHI          -8(%rbp)
#define  NMID         -16(%rbp)
#define  N            -24(%rbp)

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
         mov          %rsi, N
         mov          N, %rdi
         call         get_forward_twiddles
         mov          %rax, C
         mov          %rdx, S
         mov          $1, N2
         mov          N2, NHI
         mov          X, %rdi
         mov          N, %rsi
         call         trivial_butterfly
         mov          NHI, %rax
         shl          $2, %rax
         mov          %rax, NHI
         shl          $2, N2
         mov          N, %rcx
         bsf          N2, %rax
         bsf          NHI, %rdx
         bsf          %rcx, %rcx
         add          $4, %rax
         add          %rdx, %rax
         sub          %rax, %rcx
         xor          %rax, %rax
         mov          $1, %rdx
         shl          %cl, %rdx
         cmp          $0, %rcx
         cmovl        %rax, %rdx
         mov          %rdx, NMID
L0:
         cmp          $0, NMID
         je           L2
         mov          X, %rdi
         mov          C, %rsi
         mov          S, %rdx
         mov          NHI, %rcx
         mov          NMID, %r8
         mov          N2, %r9
         call         butterfly4_and_transpose
         mov          NHI, %rax
         shl          $2, %rax
         mov          %rax, NHI
         shl          $2, N2
         mov          N, %rcx
         bsf          N2, %rax
         bsf          NHI, %rdx
         bsf          %rcx, %rcx
         add          $4, %rax
         add          %rdx, %rax
         sub          %rax, %rcx
         xor          %rax, %rax
         mov          $1, %rdx
         shl          %cl, %rdx
         cmp          $0, %rcx
         cmovl        %rax, %rdx
         mov          %rdx, NMID
         jmp          L0
L2:
         mov          N2, %rax
         mov          N, %rcx
         shl          $2, %rax
         shr          $2, %rcx
         cmp          %rax, %rcx
         jl           L3
         mov          N, %rcx
         bsf          %rcx, %rcx
         bsf          N2, %rax
         add          $2, %rax
         sub          %rax, %rcx
         mov          $1, %rdx
         shl          %cl, %rdx
         mov          %rdx, NHI
         mov          X, %rdi
         mov          C, %rsi
         mov          S, %rdx
         mov          NHI, %rcx
         mov          N2, %r8
	     call         butterfly4
         shl          $2, N2
         jmp          L2
L3:
         mov          X, %rdi
         mov          N, %rsi
         shr          $4, %rsi
         call         transpose
         mov          X, %rdi
         mov          C, %rsi
         mov          S, %rdx
         mov          N, %rcx
         mov          N2, %r8
         call         butterfly4_finish
         pop          %r15
         pop          %r14
         pop          %r13
         pop          %r12
         pop          %rbx
         leave
         ret
