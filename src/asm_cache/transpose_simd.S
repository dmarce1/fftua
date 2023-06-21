#define  X            %rdi
#define  N            %rsi
#define  Nm1          %r8
#define  No2          %r9
#define  I            %rbx
#define  J            %rcx
#define  K            %rdx

         .global      fft_transpose_simd

         .text

fft_transpose_simd:
         push         %rbx
         shr          $2, N
         mov          N, Nm1
         dec          Nm1
         mov          N, No2
         shr          No2
         xor          I, I
         xor          J, J
L0:      cmp          I, J
         jle          L5
         vmovapd      (X, I, 8), %ymm0
         vmovapd      (X, J, 8), %ymm1
         vmovapd      %ymm1, (X, J, 8)
         vmovapd      %ymm0, (X, I, 8)
L5:      mov          No2, K
L10:     cmp          K, J
         jl           L20
         sub          K, J
         shr          K
         jmp          L10
L20:     add          K, J
         cmp          J, Nm1
         jne          L0

         pop          %rbx
         ret
