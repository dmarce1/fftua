#define  N1           $4
#define  X            %rdi
#define  N2           %rsi
#define  K            %r8
#define  n1           %r9
#define  n2           %r10
#define  I            %r11
#define  J            %r12
#define  N1N2         %r13

         .global      transpose

         .text

transpose:
         push         %r13
         push         %r12
         push         %rbx
         mov          N2, N1N2
         shl          $2, N1N2
         xor          K, K
loop:
         mov          K, %rsi
         mov          K, %rdi
         inc          %rdi
         add          N1N2, %rsi
         movq         (X, %rdi, 8), %xmm0
         movq         (X, %rsi, 8), %xmm6
         inc          %rdi
         add          N1N2, %rsi
         movq         (X, %rdi, 8), %xmm1
         movq         (X, %rsi, 8), %xmm7
         inc          %rdi
         add          N1N2, %rsi
         movq         (X, %rdi, 8), %xmm2
         movq         (X, %rsi, 8), %xmm8
         dec          %rdi
         sub          N1N2, %rsi
         add          N1N2, %rdi
         inc          %rsi
         movq         (X, %rdi, 8), %xmm3
         movq         (X, %rsi, 8), %xmm9
         inc          %rdi
         add          N1N2, %rsi
         movq         (X, %rdi, 8), %xmm4
         movq         (X, %rsi, 8), %xmm10
         add          N1N2, %rdi
         inc          %rsi
         movq         (X, %rdi, 8), %xmm5
         movq         (X, %rsi, 8), %xmm11
         mov          K, %rsi
         mov          K, %rdi
         inc          %rdi
         add          N1N2, %rsi
         movq         %xmm6, (X, %rdi, 8)
         movq         %xmm0, (X, %rsi, 8)
         inc          %rdi
         add          N1N2, %rsi
         movq         %xmm7, (X, %rdi, 8)
         movq         %xmm1, (X, %rsi, 8)
         inc          %rdi
         add          N1N2, %rsi
         movq         %xmm8, (X, %rdi, 8)
         movq         %xmm2, (X, %rsi, 8)
         dec          %rdi
         sub          N1N2, %rsi
         add          N1N2, %rdi
         inc          %rsi
         movq         %xmm9, (X, %rdi, 8)
         movq         %xmm3, (X, %rsi, 8)
         inc          %rdi
         add          N1N2, %rsi
         movq         %xmm10, (X, %rdi, 8)
         movq         %xmm4, (X, %rsi, 8)
         add          N1N2, %rdi
         inc          %rsi
         movq         %xmm11, (X, %rdi, 8)
         movq         %xmm5, (X, %rsi, 8)
         add          N1, K
         cmp          N1N2, K
         jne          loop
         pop          %rbx
         pop          %r12
         pop          %r13
         ret

