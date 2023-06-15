#define  I            %r8
#define  J            %r9
#define  K            %r10
#define  X            %r11
#define  A            %r12
#define  N1N2         %r14
#define  N1           %r15
#define  Mask         %ymm9
#define  Indices      %ymm8
#define  N2           -8(%rbp)
#define  Aptr         -1032(%rbp)
#define  STACK_SIZE   $1046

         .global      fft_transpose_hilo

         .data

         .align       32
indices: .quad        0
         .quad        4
         .quad        8
         .quad        12
mask:    .double      -1.0
         .double      -1.0
         .double      -1.0
         .double      -1.0

         .text

fft_transpose_hilo:
         enter        STACK_SIZE, $0
         push         %rbx
         push         %r12
         push         %r13
         push         %r14
         push         %r15
         mov          %rdi, X
         mov          %rsi, N1
         mov          %rdx, N2
         imul         %rsi, %rdx
         mov          %rdx, N1N2
         lea          Aptr, A
         and          $0xffffffffffffffe0, A
         vmovapd      mask, Mask
         vmovdqa      indices, Indices
         xor          I, I
L0:      xor          J, J
L10:     mov          I, %rsi
         mov          J, %rax
         imul         N1N2, %rsi
         imul         N1, %rax
         add          %rax, %rsi
         add          I, %rsi
         lea          (X, %rsi, 8), %rax
         lea          (%rax, N1N2, 8), %rbx
         lea          (%rbx, N1N2, 8), %rcx
         lea          (%rcx, N1N2, 8), %rdx
         vmovapd      (%rax), %ymm0
         vmovapd      (%rbx), %ymm1
         vmovapd      (%rcx), %ymm2
         vmovapd      (%rdx), %ymm3
         vmovapd      %ymm0, (A)
         vmovapd      %ymm1, 32(A)
         vmovapd      %ymm2, 64(A)
         vmovapd      %ymm3, 96(A)
         vmovapd      Mask, %ymm10
         vmovapd      Mask, %ymm11
         vmovapd      Mask, %ymm12
         vmovapd      Mask, %ymm13
         vgatherqpd   %ymm10, (A, Indices, 8), %ymm0
         vgatherqpd   %ymm11, 8(A, Indices, 8), %ymm1
         vgatherqpd   %ymm12, 16(A, Indices, 8), %ymm2
         vgatherqpd   %ymm13, 24(A, Indices, 8), %ymm3
         vmovapd      %ymm0, (%rax)
         vmovapd      %ymm1, (%rbx)
         vmovapd      %ymm2, (%rcx)
         vmovapd      %ymm3, (%rdx)
         mov          I, K
L20:     add          $4, K
         cmp          K, N1
         jle          L100
         mov          I, %rsi
         mov          J, %rax
         imul         N1N2, %rsi
         imul         N1, %rax
         add          %rax, %rsi
         add          K, %rsi
         mov          K, %rdi
         mov          J, %rax
         imul         N1N2, %rdi
         imul         N1, %rax
         add          %rax, %rdi
         add          I, %rdi
         lea          (X, %rsi, 8), %rax
         lea          (%rax, N1N2, 8), %rbx
         lea          (%rbx, N1N2, 8), %rcx
         lea          (%rcx, N1N2, 8), %rdx
         vmovapd      (%rax), %ymm0
         vmovapd      (%rbx), %ymm1
         vmovapd      (%rcx), %ymm2
         vmovapd      (%rdx), %ymm3
         push         %rdx
         push         %rcx
         push         %rbx
         push         %rax
         lea          (X, %rdi, 8), %rax
         lea          (%rax, N1N2, 8), %rbx
         lea          (%rbx, N1N2, 8), %rcx
         lea          (%rcx, N1N2, 8), %rdx
         vmovapd      (%rax), %ymm4
         vmovapd      (%rbx), %ymm5
         vmovapd      (%rcx), %ymm6
         vmovapd      (%rdx), %ymm7
         vmovapd      %ymm0, (A)
         vmovapd      %ymm1, 32(A)
         vmovapd      %ymm2, 64(A)
         vmovapd      %ymm3, 96(A)
         vmovapd      %ymm4, 128(A)
         vmovapd      %ymm5, 160(A)
         vmovapd      %ymm6, 192(A)
         vmovapd      %ymm7, 224(A)
         vmovapd      Mask, %ymm10
         vmovapd      Mask, %ymm11
         vmovapd      Mask, %ymm12
         vmovapd      Mask, %ymm13
         vgatherqpd   %ymm10, (A, Indices, 8), %ymm0
         vgatherqpd   %ymm11, 8(A, Indices, 8), %ymm1
         vgatherqpd   %ymm12, 16(A, Indices, 8), %ymm2
         vgatherqpd   %ymm13, 24(A, Indices, 8), %ymm3
         vmovapd      Mask, %ymm10
         vmovapd      Mask, %ymm11
         vmovapd      Mask, %ymm12
         vmovapd      Mask, %ymm13
         vgatherqpd   %ymm10, 128(A, Indices, 8), %ymm4
         vgatherqpd   %ymm11, 136(A, Indices, 8), %ymm5
         vgatherqpd   %ymm12, 144(A, Indices, 8), %ymm6
         vgatherqpd   %ymm13, 152(A, Indices, 8), %ymm7
         vmovntpd     %ymm0, (%rax)
         vmovntpd     %ymm1, (%rbx)
         vmovntpd     %ymm2, (%rcx)
         vmovntpd     %ymm3, (%rdx)
         pop          %rax
         pop          %rbx
         pop          %rcx
         pop          %rdx
         vmovntpd     %ymm4, (%rax)
         vmovntpd     %ymm5, (%rbx)
         vmovntpd     %ymm6, (%rcx)
         vmovntpd     %ymm7, (%rdx)
         jmp          L20
L100:    inc          J
         cmp          J, N2
         jne          L10
         add          $4, I
         cmp          I, N1
         jne          L0
         pop          %r15
         pop          %r14
         pop          %r13
         pop          %r12
         pop          %rbx
         leave
         ret

