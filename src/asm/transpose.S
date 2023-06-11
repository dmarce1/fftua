#define  N1           $4
#define  X            %rdi
#define  N2           %rsi
#define  n0           %r8
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
         xor          n0, n0
loop0:
         xor          n1, n1
loop1:
         mulx         n0, N2, I
         add          n1, I
         shl          $5, I
         mov          n1, J
         shl          $2, J
         add          n0, J
         lea          (X, J, 8), %rax
         lea          (%rax, N1N2, 8), %rbx
         lea          (%rbx, N1N2, 8), %rcx
         lea          (%rcx, N1N2, 8), %rdx
         movq         0(X, I), %xmm0
         movq         8(X, I), %xmm1
         movq         16(X, I), %xmm2
         movq         24(X, I), %xmm3
         movq         (%rax), %xmm4
         movq         (%rbx), %xmm5
         movq         (%rcx), %xmm6
         movq         (%rdx), %xmm7
         movq         %xmm4, 0(X, I)
         movq         %xmm5, 8(X, I)
         movq         %xmm6, 16(X, I)
         movq         %xmm7, 24(X, I)
         movq         %xmm0, (%rax)
         movq         %xmm1, (%rbx)
         movq         %xmm2, (%rcx)
         movq         %xmm3, (%rdx)
         inc          n1
         cmp          N2, n1
         jne          loop1
         inc          n0
         cmp          N1, n0
         jne          loop0
         pop          %rbx
         pop          %r12
         pop          %r13
         ret

