#define  NHIo2        %rbx
#define  NHIm1        %rcx
#define  I0           %rsi
#define  J0           %rdi
#define  I            %r8
#define  J            %r9
#define  K            %r10
#define  L            %r11
#define  X            %r12
#define  NHI          %r13
#define  NLO          %r15

         .global      fft_scramble_hi

         .text

fft_scramble_hi:
         push         %rbx
         push         %r12
         push         %r13
         push         %r14
         push         %r15
         mov          %rdi, X
         mov          %rsi, NHI
         mov          %rdx, NLO
         cmp          $2, NHI
         jle          skip
         mov          NHI, NHIm1
         mov          NHI, NHIo2
         shr          NHIo2
         dec          NHIm1
         xor          I, I
         xor          J, J
L0:      cmp          I, J
         jg           L5
         mov          I, I0
         mov          J, J0
         imul         NLO, I0
         imul         NLO, J0
         mov          I0, %rax
         mov          J0, %rdx
         mov          I0, L
         add          NLO, L
L1:      vmovapd      (X, %rax, 8), %ymm0
         vmovapd      (X, %rdx, 8), %ymm1
         vmovapd      %ymm1, (X, %rax, 8)
         vmovapd      %ymm0, (X, %rdx, 8)
         add          $4, %rax
         add          $4, %rdx
         cmp          %rax, L
         jg           L1
L5:      mov          NHIo2, K
L10:     cmp          K, J
         jl           L20
         sub          K, J
         shr          K
         jmp          L10
L20:     add          K, J
         inc          I
         cmp          I, NHIm1
         jg           L0
skip:    pop          %r15
         pop          %r14
         pop          %r13
         pop          %r12
         pop          %rbx
         ret

