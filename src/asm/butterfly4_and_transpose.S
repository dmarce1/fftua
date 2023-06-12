#define  STACK_SIZE   $1216
#define  M_SQRT1_2    0.70710678118654752440
#define  N1           $4
#define  X            %r8
#define  W            %r9
#define  imid         %r11
#define  k2           %r12
#define  n0           %r13
#define  n1           %r14
#define  U0           %r15
#define  er0          %ymm0
#define  er1          %ymm1
#define  er2          %ymm2
#define  er3          %ymm3
#define  ei0          %ymm4
#define  ei1          %ymm5
#define  ei2          %ymm6
#define  ei3          %ymm7
#define  tr0          %ymm8
#define  tr1          %ymm9
#define  tr2          %ymm10
#define  tr3          %ymm11
#define  ti0          %ymm12
#define  ti1          %ymm13
#define  ti2          %ymm14
#define  ti3          %ymm15
#define  NHI          -8(%rbp)
#define  N2           -16(%rbp)
#define  NMID         -24(%rbp)
#define  NHIoN1       -32(%rbp)
#define  N1NMID       -40(%rbp)
#define  N1N2NMID     -48(%rbp)
#define  N            -56(%rbp)
#define  TWHI         -64(%rbp)
#define  N2o2         -72(%rbp)
#define  cos1         -80(%rbp)
#define  cos2         -88(%rbp)
#define  cos3         -96(%rbp)
#define  sin1         -104(%rbp)
#define  sin2         -112(%rbp)
#define  sin3         -120(%rbp)
#define  ihi          -128(%rbp)
#define  klor         -136(%rbp)
#define  khir         -144(%rbp)
#define  kloi         -152(%rbp)
#define  khii         -160(%rbp)
#define  U0ptr        -1184(%rbp)

         .global      butterfly4_and_transpose

         .data

         .align       32
tw45:    .double      M_SQRT1_2
         .double      M_SQRT1_2
         .double      M_SQRT1_2
         .double      M_SQRT1_2
none:    .double      -1.0
         .double      -1.0
         .double      -1.0
         .double      -1.0
two:     .quad        2
three:   .quad        3

         .text

butterfly4_and_transpose:
         enter        STACK_SIZE, $0
         push         %rbx
         push         %r12
         push         %r13
         push         %r14
         push         %r15
         mov          %r8, N2
         mov          %rcx, NMID
         mov          %rdx, NHI
         mov          %rsi, W
         mov          %rdi, X
         lea          U0ptr, U0
         and          0xffffffffffffffe0, U0
         mov          NHI, %rax
         shr          $2, %rax
         mov          %rax, NHIoN1
         mov          NMID, %rbx
         shl          $2, %rbx
         mov          %rbx, N1NMID
         mov          N2, %rax
         mul          %rbx
         mov          %rax, N1N2NMID
         mov          NHI, %rbx
         mul          %rbx
         shl          $2, %rax
         mov          %rax, N
         mov          N, %rax
         mov          N2, %rbx
         shl          $2, %rbx
         div          %rbx
         mov          %rax, TWHI
         xor          %rdx, %rdx
         mov          %rdx, ihi
ihiloop:
         xor          imid, imid
imidloop0:
         mov          N1NMID, %rax
         mov          ihi, %rbx
         mul          %rbx
         add          imid, %rax
         mov          N2, %rbx
         mul          %rbx
         shl          $2, %rax
         mov          %rax, %rdi
         xor          n0, n0
n0loop2:
         xor          n1, n1
n1loop2:
         mov          n1, %rax
         mov          N1NMID, %rbx
         mul          %rbx
         add          n0, %rax
         mov          N2, %rbx
         mul          %rbx
         mov          %rax, %rdx
         add          %rdx, %rdi
         mov          n0, %rax
         shl          $2, %rax
         add          n1, %rax
         shl          $2, %rax
         vmovapd      (X, %rdx, 8), %ymm0
         vmovupd      %ymm0, (U0, %rax, 8)
         inc          n1
         cmp          N1, n1
         jne          n1loop2
         inc          n0
         cmp          N1, n0
         jne          n0loop2
         xor          n0, n0
n0loop3:
         mov          n0, %rdx
         shl          $7, %rdx
         add          U0, %rdx
         vmovupd      (%rdx), er0
         vmovupd      32(%rdx), er1
         vmovupd      64(%rdx), er2
         vmovupd      96(%rdx), er3
         vaddpd       er2, er0, tr0
         vsubpd       er2, er0, tr2
         vaddpd       er3, er1, tr1
         vsubpd       er3, er1, tr3
         vaddpd       tr1, tr0, er0
         vsubpd       tr1, tr0, er2
         mov          n0, %rax
         mov          N1N2NMID, %rbx
         mul          %rbx
         add          %rdi, %rax
         mov          N2, %rdx
         lea          (X, %rax, 8), %rax
         vmovapd      er0, (%rax)
         lea          (%rax, %rdx, 8), %rax
         vmovapd      tr2, (%rax)
         lea          (%rax, %rdx, 8), %rax
         vmovapd      tr3, (%rax)
         lea          (%rax, %rdx, 8), %rax
         vmovapd      er2, (%rax)
         inc          n0
         cmp          N1, n0
         jne          n0loop3
         inc          imid
         cmp          imid, NMID
         jne          imidloop0
         xor          imid, imid
imidloop1:
         mov          N2o2, %rax
         mov          N2o2, %rbx
         shr          $2, %rax
         and          $3, %rbx
         mov          %rax, khir
         mov          %rbx, klor
         mov          NHIoN1, %rax
         mov          klor, %rbx
         mul          %rbx
         add          ihi, %rax
         mov          N1NMID, %rbx
         mul          %rbx
         add          imid, %rax
         mov          N2, %rbx
         mul          %rbx
         add          khir, %rax
         shl          $2, %rax
         mov          %rax, %rdi
         xor          n0, n0
n0loop4:
         xor          n1, n1
n1loop4:
         mov          n1, %rax
         mov          N1NMID, %rbx
         mul          %rbx
         add          n0, %rax
         mov          N2, %rbx
         mul          %rbx
         mov          %rax, %rdx
         add          %rdi, %rdx
         mov          n0, %rax
         shl          $2, %rax
         add          n1, %rax
         shl          $2, %rax
         vmovapd      (X, %rdx, 8), %ymm0
         vmovupd      %ymm0, (U0, %rax, 8)
         inc          n1
         cmp          N1, n1
         jne          n1loop2
         inc          n0
         cmp          N1, n0
         jne          n0loop2
         xor          n0, n0
n0loop5:
         mov          n0, %rdx
         shl          $7, %rdx
         add          U0, %rdx
         vmovupd      (%rdx), er0
         vmovupd      32(%rdx), er1
         vmovupd      64(%rdx), er2
         vmovupd      96(%rdx), er3
         vaddpd       er3, er1, tr0
         vsubpd       er3, er1, tr2
         vmulpd       tw45, tr2, tr1
         vmulpd       tw45, tr0, tr3
         vmovapd      er0, tr0
         vmovapd      er2, tr2
         vaddpd       tr1, tr0, er0
         vaddpd       tr3, tr2, er3
         vsubpd       tr1, tr0, er1
         vsubpd       tr3, tr2, er2
         vmulpd       none, er3, er3
         mov          N1N2NMID, %rax
         mov          n0, %rbx
         mul          %rax
         add          %rdi, %rax
         mov          N2, %rdx
         lea          (X, %rax, 8), %rax
         vmovapd      er0, (%rax)
         lea          (%rax, %rdx, 8), %rax
         vmovapd      er1, (%rax)
         lea          (%rax, %rdx, 8), %rax
         vmovapd      er2, (%rax)
         lea          (%rax, %rdx, 8), %rax
         vmovapd      er3, (%rax)
         inc          n0
         cmp          N1, n0
         jne          n0loop5
         inc          imid
         cmp          imid, NMID
         jne          imidloop1
imidloop2:
         xor          k2, k2
k2loop_b:
         inc          k2
         cmp          k2, N2o2
         je           k2loop_e
         mov          k2, %rbx
         mov          TWHI, %rax
         mul          %rbx
         shl          %rax
         mov          %rax, %rbx
         mov          %rax, %rcx
         add          %rax, %rbx
         add          %rbx, %rcx
         movq         0(W, %rax, 8), %xmm0
         movq         8(W, %rax, 8), %xmm1
         movq         0(W, %rbx, 8), %xmm2
         movq         8(W, %rbx, 8), %xmm3
         movq         0(W, %rcx, 8), %xmm4
         movq         8(W, %rcx, 8), %xmm5
         movq         %xmm0, cos1
         movq         %xmm1, sin1
         movq         %xmm2, cos2
         movq         %xmm3, sin2
         movq         %xmm4, cos3
         movq         %xmm5, sin3
         mov          k2, %rax
         mov          k2, %rbx
         mov          N2, %rcx
         sub          k2, %rcx
         mov          %rcx, %rdx
         shr          $2, %rax
         and          $3, %rbx
         shr          $2, %rcx
         and          $3, %rdx
         mov          %rax, khir
         mov          %rbx, klor
         mov          %rcx, khii
         mov          %rdx, kloi
         mov          NHIoN1, %rax
         mov          klor, %rbx
         mul          %rbx
         add          ihi, %rax
         mov          N1NMID, %rbx
         mul          %rbx
         add          imid, %rax
         mov          N2, %rbx
         mul          %rbx
         add          khir, %rax
         mov          %rax, %rdi
         mov          NHIoN1, %rax
         mov          kloi, %rbx
         mul          %rbx
         add          ihi, %rax
         mov          N1NMID, %rbx
         mul          %rbx
         add          imid, %rax
         mov          N2, %rbx
         mul          %rbx
         add          khii, %rax
         mov          %rax, %rsi
         shl          $2, %rdi
         shl          $2, %rsi
         xor          n0, n0
n0loop0:
         xor          n1, n1
n1loop0:
         mov          n1, %rax
         mov          N1NMID, %rbx
         mul          %rbx
         add          n0, %rax
         mov          N2, %rbx
         mul          %rbx
         mov          %rax, %rcx
         mov          %rax, %rdx
         add          %rdi, %rcx
         add          %rsi, %rdx
         vmovdqa      (X, %rcx, 8), %ymm0
         vmovdqa      (X, %rdx, 8), %ymm1
         mov          n0, %rax
         shl          $2, %rax
         add          n1, %rax
         shl          $3, %rax
         vmovdqu      %ymm0, 0(U0, %rax, 8)
         vmovdqu      %ymm1, 32(U0, %rax, 8)
         inc          n1
         cmp          N1, n1
         jne          n1loop0
         inc          n0
         cmp          N1, n0
         jne          n0loop0
         xor          n0, n0
n0loop1:
         mov          n0, %rdx
         shl          $8, %rdx
         add          U0, %rdx
         vbroadcastsd cos1, tr1
         vbroadcastsd cos2, tr2
         vbroadcastsd cos3, tr3
         vbroadcastsd sin1, ti1
         vbroadcastsd sin2, ti2
         vbroadcastsd sin3, ti3
         vmovupd      (%rdx), er0
         vmovupd      32(%rdx), ei0
         vmovupd      64(%rdx), er1
         vmovupd      96(%rdx), ei1
         vmovupd      128(%rdx), er2
         vmovupd      160(%rdx), ei2
         vmovupd      192(%rdx), er3
         vmovupd      224(%rdx), ei3
         vmulpd       ti1, ei1, tr0
         vmulpd       tr1, ei1, ti0
         vfmsub231pd  tr1, er1, tr0
         vfmadd231pd  ti1, er1, ti0
         vmovapd      ti0, ei1
         vmovapd      tr0, er1
         vmulpd       ti2, ei2, tr0
         vmulpd       tr2, ei2, ti0
         vfmsub231pd  tr2, er2, tr0
         vfmadd231pd  ti2, er2, ti0
         vmovapd      ti0, ei2
         vmovapd      tr0, er2
         vmulpd       ti3, ei3, tr0
         vmulpd       tr3, ei3, ti0
         vfmsub231pd  tr3, er3, tr0
         vfmadd231pd  ti3, er3, ti0
         vmovapd      ti0, ei3
         vmovapd      tr0, er3
         vaddpd       er2, er0, tr0
         vaddpd       ei2, ei0, ti0
         vsubpd       er2, er0, tr2
         vsubpd       ei2, ei0, ti2
         vaddpd       er3, er1, tr1
         vaddpd       ei3, ei1, ti1
         vsubpd       er1, er3, tr3
         vsubpd       ei3, ei1, ti3
         vaddpd       tr1, tr0, er0
         vaddpd       ti1, ti0, ei0
         vaddpd       ti3, tr2, er1
         vaddpd       tr3, ti2, ei1
         vsubpd       tr1, tr0, er2
         vsubpd       ti0, ti1, ei2
         vsubpd       ti3, tr2, er3
         vsubpd       ti2, tr3, ei3
         mov          n0, %rax
         mov          N1N2NMID, %rbx
         mul          %rbx
         add          %rax, %rdi
         add          %rax, %rsi
         mov          N2, %rdx
         lea          (X, %rdi, 8), %rax
         lea          (X, %rsi, 8), %rbx
         vmovapd      er0, (%rax)
         vmovapd      er1, (%rax)
         lea          (%rax, %rdx, 8), %rax
         lea          (%rbx, %rdx, 8), %rbx
         vmovapd      ei2, (%rax)
         vmovapd      ei3, (%rbx)
         lea          (%rax, %rdx, 8), %rax
         lea          (%rbx, %rdx, 8), %rbx
         vmovapd      er3, (%rax)
         vmovapd      er2, (%rbx)
         lea          (%rax, %rdx, 8), %rax
         lea          (%rbx, %rdx, 8), %rbx
         vmovapd      ei1, (%rax)
         vmovapd      ei0, (%rbx)
         inc          n0
         cmp          N1, n0
         jne          n0loop1
         jmp          k2loop_b
k2loop_e:
         inc          imid
         cmp          imid, NMID
         jne          imidloop2
         mov          ihi, %rax
         inc          %rax
         mov          %rax, ihi
         cmp          %rax, NHIoN1
         jne          ihiloop
         pop          %r15
         pop          %r14
         pop          %r13
         pop          %r12
         pop          %rbx
         leave
         ret

