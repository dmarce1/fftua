#define  I            %r8
#define  J            %r9
#define  K            %r10
#define  X            %r11
#define  A            %r12
#define  N1N2         %r14
#define  N1           %r15
#define  Mask         %ymm14
#define  Indices      %ymm15
#define  N2           -8(%rbp)
#define  tmp3         -16(%rbp)
#define  tmp2         -24(%rbp)
#define  tmp1         -32(%rbp)
#define  tmp0         -40(%rbp)
#define  Aptr         -552(%rbp)
#define  STACK_SIZE   $576

         .global      fft_transpose_hilo

         .data

         .align       32
mask:    .double      -1.0
         .double      -1.0
         .double      -1.0
         .double      -1.0
indices: .long        0
         .long        4
         .long        8
         .long        12

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
         xor          %rax, %rax
         mov          %rax, tmp0
         add          %rdx, %rax
         mov          %rax, tmp1
         add          %rdx, %rax
         mov          %rax, tmp2
         add          %rdx, %rax
         mov          %rax, tmp3
         vmovdqu      tmp0, Indices
         lea          Aptr, A
         and          $0xffffffffffffffe0, A
         vmovapd      mask, Mask
         xor          I, I
L0:      xor          J, J
L10:     mov          I, %rsi
         mov          J, %rax
         imul         N1N2, %rsi
         imul         N1, %rax
         add          %rax, %rsi
         add          I, %rsi
         lea          (X, %rsi, 8), %rsi
         vmovapd      Mask, %ymm4
         vmovapd      Mask, %ymm5
         vmovapd      Mask, %ymm6
         vmovapd      Mask, %ymm7
         vgatherqpd   %ymm4, (%rsi, Indices, 8), %ymm0
         vgatherqpd   %ymm5, 8(%rsi, Indices, 8), %ymm1
         vgatherqpd   %ymm6, 16(%rsi, Indices, 8), %ymm2
         vgatherqpd   %ymm7, 24(%rsi, Indices, 8), %ymm3
         lea          (%rsi), %rax
         lea          (%rax, N1N2, 8), %rbx
         lea          (%rbx, N1N2, 8), %rcx
         lea          (%rcx, N1N2, 8), %rdx
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
         lea          (X, %rsi, 8), %rsi
         lea          (X, %rdi, 8), %rdi
         vmovapd      Mask, %ymm4
         vmovapd      Mask, %ymm5
         vmovapd      Mask, %ymm6
         vmovapd      Mask, %ymm7
         vmovapd      Mask, %ymm8
         vmovapd      Mask, %ymm9
         vmovapd      Mask, %ymm10
         vmovapd      Mask, %ymm11
         vgatherqpd   %ymm4, (%rdi, Indices, 8), %ymm0
         vgatherqpd   %ymm5, 8(%rdi, Indices, 8), %ymm1
         vgatherqpd   %ymm6, 16(%rdi, Indices, 8), %ymm2
         vgatherqpd   %ymm7, 24(%rdi, Indices, 8), %ymm3
         vgatherqpd   %ymm8, (%rsi, Indices, 8), %ymm4
         vgatherqpd   %ymm9, 8(%rsi, Indices, 8), %ymm5
         vgatherqpd   %ymm10, 16(%rsi, Indices, 8), %ymm6
         vgatherqpd   %ymm11, 24(%rsi, Indices, 8), %ymm7
         lea          (%rsi), %rax
         lea          (%rax, N1N2, 8), %rbx
         lea          (%rbx, N1N2, 8), %rcx
         lea          (%rcx, N1N2, 8), %rdx
         vmovapd      %ymm0, (%rax)
         vmovapd      %ymm1, (%rbx)
         vmovapd      %ymm2, (%rcx)
         vmovapd      %ymm3, (%rdx)
         lea          (%rdi), %rax
         lea          (%rax, N1N2, 8), %rbx
         lea          (%rbx, N1N2, 8), %rcx
         lea          (%rcx, N1N2, 8), %rdx
         vmovapd      %ymm4, (%rax)
         vmovapd      %ymm5, (%rbx)
         vmovapd      %ymm6, (%rcx)
         vmovapd      %ymm7, (%rdx)
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






