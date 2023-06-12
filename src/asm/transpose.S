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
         movq         (X, %rdx, 8), %ymm0
         movq         (X, %rax, 8), %ymm6
         inc          %rdx
         add          N1N2, %rax
         movq         (X, %rdx, 8), %ymm1
         movq         (X, %rax, 8), %ymm7
         inc          %rdx
         add          N1N2, %rax
         movq         (X, %rdx, 8), %ymm2
         movq         (X, %rax, 8), %ymm8
         dec          %rdx
         sub          N1N2, %rax
         add          N1N2, %rdx
         inc          %rax
         movq         (X, %rdx, 8), %ymm3
         movq         (X, %rax, 8), %ymm9
         inc          %rdx
         add          N1N2, %rax
         movq         (X, %rdx, 8), %ymm4
         movq         (X, %rax, 8), %ymm10
         add          N1N2, %rdx
         inc          %rax
         movq         (X, %rdx, 8), %ymm5
         movq         (X, %rax, 8), %ymm11

         mov          K, %rax
         mov          K, %rdx
         inc          %rdx
         add          N1N2, %rax
         movq         %ymm6, (X, %rdx, 8)
         movq         %ymm0, (X, %rax, 8)
         inc          %rdx
         add          N1N2, %rax
         movq         %ymm7, (X, %rdx, 8)
         movq         %ymm1, (X, %rax, 8)
         inc          %rdx
         add          N1N2, %rax
         movq         %ymm8, (X, %rdx, 8)
         movq         %ymm2, (X, %rax, 8)
         dec          %rdx
         sub          N1N2, %rax
         add          N1N2, %rdx
         inc          %rax
         movq         %ymm9, (X, %rdx, 8)
         movq         %ymm3, (X, %rax, 8)
         inc          %rdx
         add          N1N2, %rax
         movq         %ymm10, (X, %rdx, 8)
         movq         %ymm4, (X, %rax, 8)
         add          N1N2, %rdx
         inc          %rax
         movq         %ymm11, (X, %rdx, 8)
         movq         %ymm5, (X, %rax, 8)
         add          N1, K
         cmp          N1N2, K
         jne          loop
         pop          %r12
         ret

