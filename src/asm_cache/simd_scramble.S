#define  I            %r8
#define  J            %r9
#define  K            %r10
#define  X            %r11
#define  N2           %r12
#define  N2o2         %r13
#define  N2m1         %r14

         .global      fft_simd_scramble

         .text

fft_simd_scramble:
         push         %r12
         push         %r13
         push         %r14
         mov          %rdi, X
         mov          %rsi, N2
         shr          $2, N2
         mov          N2, N2m1
         mov          N2, N2o2
         shr          N2o2
         dec          N2m1
         xor          I, I
         xor          J, J
L0:      cmp          I, J
         jg           L5
         mov          I, %rax
         mov          J, %rdx
         shl          $5, %rax
         shl          $5, %rdx
         vmovapd      (X, %rax), %ymm0
         vmovapd      (X, %rdx), %ymm1
         vmovapd      %ymm0, (X, %rdx)
         vmovapd      %ymm1, (X, %rax)
L5:      mov          N2o2, K
L10:     cmp          K, J
         jl           L20
         sub          K, J
         shr          K
         jmp          L10
L20:     add          K, J
         inc          I
         cmp          I, N2m1
         jg           L0
         pop          %r14
         pop          %r13
         pop          %r12
         ret

