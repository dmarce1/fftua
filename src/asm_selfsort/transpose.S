#define  N1           $4
#define  X            %rdi
#define  N1N2         %rsi
#define  K            %r8

         .global      transpose

         .text

transpose:
         shl          $2, N1N2
         xor          K, K
loop:
         mov          K, %rax
         mov          K, %rdx
         inc          %rdx
         add          N1N2, %rax
         mov          (X, %rdx, 8), %r10
         mov          (X, %rax, 8), %r11
         mov          %r11, (X, %rdx, 8)
         mov          %r10, (X, %rax, 8)
         inc          %rdx
         add          N1N2, %rax
         mov          (X, %rdx, 8), %r10
         mov          (X, %rax, 8), %r11
         mov          %r11, (X, %rdx, 8)
         mov          %r10, (X, %rax, 8)
         inc          %rdx
         add          N1N2, %rax
         mov          (X, %rdx, 8), %r10
         mov          (X, %rax, 8), %r11
         mov          %r11, (X, %rdx, 8)
         mov          %r10, (X, %rax, 8)
         dec          %rdx
         sub          N1N2, %rax
         add          N1N2, %rdx
         inc          %rax
         mov          (X, %rdx, 8), %r10
         mov          (X, %rax, 8), %r11
         mov          %r11, (X, %rdx, 8)
         mov          %r10, (X, %rax, 8)
         inc          %rdx
         add          N1N2, %rax
         mov          (X, %rdx, 8), %r10
         mov          (X, %rax, 8), %r11
         mov          %r11, (X, %rdx, 8)
         mov          %r10, (X, %rax, 8)
         add          N1N2, %rdx
         inc          %rax
         mov          (X, %rdx, 8), %r10
         mov          (X, %rax, 8), %r11
         mov          %r11, (X, %rdx, 8)
         mov          %r10, (X, %rax, 8)
         add          N1, K
         cmp          N1N2, K
         jne          loop
         ret

