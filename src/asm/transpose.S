#define  N1           $4
#define  X            %rdi
#define  N2           %rsi
#define  K            %r8
#define  n1           %r9
#define  n2           %r10
#define  I            %r11
#define  J            %r12
#define  N1N2         %rcx

         .global      transpose

         .text

transpose:
         push         %r12
         mov          N2, N1N2
         shl          $2, N1N2
         xor          K, K
loop:
         mov          K, %rax
         mov          K, %rdx

         inc          %rdx
         add          N1N2, %rax
         movq         (X, %rdx, 8), %xmm0
         movq         (X, %rax, 8), %xmm6
         inc          %rdx
         add          N1N2, %rax
         movq         (X, %rdx, 8), %xmm1
         movq         (X, %rax, 8), %xmm7
         inc          %rdx
         add          N1N2, %rax
         movq         (X, %rdx, 8), %xmm2
         movq         (X, %rax, 8), %xmm8
         dec          %rdx
         sub          N1N2, %rax
         add          N1N2, %rdx
         inc          %rax
         movq         (X, %rdx, 8), %xmm3
         movq         (X, %rax, 8), %xmm9
         inc          %rdx
         add          N1N2, %rax
         movq         (X, %rdx, 8), %xmm4
         movq         (X, %rax, 8), %xmm10
         add          N1N2, %rdx
         inc          %rax
         movq         (X, %rdx, 8), %xmm5
         movq         (X, %rax, 8), %xmm11

         mov          K, %rax
         mov          K, %rdx
         inc          %rdx
         add          N1N2, %rax
         movq         %xmm6, (X, %rdx, 8)
         movq         %xmm0, (X, %rax, 8)
         inc          %rdx
         add          N1N2, %rax
         movq         %xmm7, (X, %rdx, 8)
         movq         %xmm1, (X, %rax, 8)
         inc          %rdx
         add          N1N2, %rax
         movq         %xmm8, (X, %rdx, 8)
         movq         %xmm2, (X, %rax, 8)
         dec          %rdx
         sub          N1N2, %rax
         add          N1N2, %rdx
         inc          %rax
         movq         %xmm9, (X, %rdx, 8)
         movq         %xmm3, (X, %rax, 8)
         inc          %rdx
         add          N1N2, %rax
         movq         %xmm10, (X, %rdx, 8)
         movq         %xmm4, (X, %rax, 8)
         add          N1N2, %rdx
         inc          %rax
         movq         %xmm11, (X, %rdx, 8)
         movq         %xmm5, (X, %rax, 8)
         add          N1, K
         cmp          N1N2, K
         jne          loop
         pop          %r12
         ret

