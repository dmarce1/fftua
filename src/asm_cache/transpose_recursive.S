#define  X            %r8
#define  M            %r9
#define  Ix           %r10
#define  Iy           %r11
#define  IMID         %r12
#define  I1           %r13
#define  I2           %r14
#define  N1N2         %r15
#define  M0           $4
#define  STACK_SIZE   $16
#define  N1           -8(%rbp)
#define  N2           -16(%rbp)

         .global      fft_transpose_recursive

         .data

shuf:    .byte        16
         .byte        17
         .byte        18
         .byte        19
         .byte        20
         .byte        21
         .byte        22
         .byte        23
         .byte        24
         .byte        25
         .byte        26
         .byte        27
         .byte        28
         .byte        29
         .byte        30
         .byte        31
         .byte         0
         .byte         1
         .byte         2
         .byte         3
         .byte         4
         .byte         5
         .byte         6
         .byte         7
         .byte         8
         .byte         9
         .byte        10
         .byte        11
         .byte        12
         .byte        13
         .byte        14
         .byte        15
         .text


fft_transpose_recursive:
         enter        STACK_SIZE, $0
         push         %r12
         push         %r13
         push         %r14
         push         %r15
         push         %rbx
         mov          %rdi, X
         mov          %rsi, N1
         mov          %rdx, N2
         imul         N1, %rdx
         mov          %rdx, N1N2
         xor          Ix, Ix
         xor          Iy, Iy
         mov          N1, M
         call         recurse
         pop          %rbx
         pop          %r15
         pop          %r14
         pop          %r13
         pop          %r12
         leave
         ret
recurse: cmp          M0, M
         jg           further
         xor          IMID, IMID
loop0:   mov          Ix, %rsi
         mov          Iy, %rdi
         imul         N2, %rsi
         imul         N2, %rdi
         add          IMID, %rsi
         add          IMID, %rdi
         imul         N1, %rsi
         imul         N1, %rdi
         add          Iy, %rsi
         add          Ix, %rdi
         lea          (X, %rsi, 8), %rsi
         lea          (X, %rdi, 8), %rdi
         cmp          %rsi, %rdi
         jl           skip
         je           ieq
         mov          N1N2, %rax
         mov          N1N2, %rbx
         add          %rax, %rax
         add          %rax, %rbx
         vmovupd      (%rsi), %ymm0
         vmovupd      (%rsi, N1N2, 8), %ymm1
         vmovupd      (%rsi, %rax, 8), %ymm2
         vmovupd      (%rsi, %rbx, 8), %ymm3
         vmovupd      (%rdi), %ymm8
         vmovupd      (%rdi, N1N2, 8), %ymm9
         vmovupd      (%rdi, %rax, 8), %ymm10
         vmovupd      (%rdi, %rbx, 8), %ymm11
         vshufpd      $0, %ymm1, %ymm0, %ymm4
         vshufpd      $15, %ymm1, %ymm0, %ymm5
         vshufpd      $0, %ymm3, %ymm2, %ymm6
         vshufpd      $15, %ymm3, %ymm2, %ymm7
         vinsertf128  $1, %xmm6, %ymm4, %ymm0
         vinsertf128  $1, %xmm7, %ymm5, %ymm1
         vpermpd      $78, %ymm4, %ymm4
         vpermpd      $78, %ymm5, %ymm5
         vinsertf128  $0, %xmm4, %ymm6, %ymm2
         vinsertf128  $0, %xmm5, %ymm7, %ymm3
         vshufpd      $0, %ymm9, %ymm8, %ymm12
         vshufpd      $15, %ymm9, %ymm8, %ymm13
         vshufpd      $0, %ymm11, %ymm10, %ymm14
         vshufpd      $15, %ymm11, %ymm10, %ymm15
         vinsertf128  $1, %xmm14, %ymm12, %ymm8
         vinsertf128  $1, %xmm15, %ymm13, %ymm9
         vpermpd      $78, %ymm12, %ymm12
         vpermpd      $78, %ymm13, %ymm13
         vinsertf128  $0, %xmm12, %ymm14, %ymm10
         vinsertf128  $0, %xmm13, %ymm15, %ymm11
         vmovupd      %ymm0, (%rdi)
         vmovupd      %ymm1, (%rdi, N1N2, 8)
         vmovupd      %ymm2, (%rdi, %rax, 8)
         vmovupd      %ymm3, (%rdi, %rbx, 8)
         vmovupd      %ymm8, (%rsi)
         vmovupd      %ymm9, (%rsi, N1N2, 8)
         vmovupd      %ymm10, (%rsi, %rax, 8)
         vmovupd      %ymm11, (%rsi, %rbx, 8)
         jmp          skip
ieq:     mov          N1N2, %rax
         mov          N1N2, %rbx
         add          %rax, %rax
         add          %rax, %rbx
         vmovupd      (%rsi), %ymm0
         vmovupd      (%rsi, N1N2, 8), %ymm1
         vmovupd      (%rsi, %rax, 8), %ymm2
         vmovupd      (%rsi, %rbx, 8), %ymm3
         vshufpd      $0, %ymm1, %ymm0, %ymm4
         vshufpd      $15, %ymm1, %ymm0, %ymm5
         vshufpd      $0, %ymm3, %ymm2, %ymm6
         vshufpd      $15, %ymm3, %ymm2, %ymm7
         vinsertf128  $1, %xmm6, %ymm4, %ymm0
         vinsertf128  $1, %xmm7, %ymm5, %ymm1
         vpermpd      $78, %ymm4, %ymm4
         vpermpd      $78, %ymm5, %ymm5
         vinsertf128  $0, %xmm4, %ymm6, %ymm2
         vinsertf128  $0, %xmm5, %ymm7, %ymm3
         vmovupd      %ymm0, (%rsi)
         vmovupd      %ymm1, (%rsi, N1N2, 8)
         vmovupd      %ymm2, (%rsi, %rax, 8)
         vmovupd      %ymm3, (%rsi, %rbx, 8)
skip:    inc          IMID
         cmp          IMID, N2
         jg           loop0
         ret
further: cmp          Ix, Iy
         je           xyeq
         shr          M
         call         recurse
         add          M, Iy
         call         recurse
         sub          M, Iy
         add          M, Ix
         call         recurse
         add          M, Iy
         call         recurse
         sub          M, Iy
         sub          M, Ix
         shl          M
         ret
xyeq:    shr          M
         call         recurse
         add          M, Iy
         call         recurse
         add          M, Ix
         call         recurse
         sub          M, Iy
         sub          M, Ix
         shl          M
         ret





