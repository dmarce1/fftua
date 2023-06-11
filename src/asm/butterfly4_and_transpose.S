#define  STACK_SIZE   $1224
#define  M_SQRT1_2    0.70710678118654752440
#define  N1           $4
#define  X            %r8
#define  W            %r9
#define  imid         %r11
#define  k2           %r12
#define  n0           %r13
#define  n1           %r14
#define  U0           %r15
#define  NHI          -8(%ebp)
#define  N2           -16(%ebp)
#define  NMID         -24(%ebp)
#define  NHIN1        -32(%ebp)
#define  N2N1         -40(%ebp)
#define  N1NMID       -48(%ebp)
#define  N1N2NMID     -56(%ebp)
#define  N            -64(%ebp)
#define  TWHI         -72(%ebp)
#define  N2o2         -80(%ebp)
#define  cos1         -88(%ebp)
#define  cos2         -96(%ebp)
#define  cos3         -104(%ebp)
#define  sin1         -112(%ebp)
#define  sin2         -120(%ebp)
#define  sin3         -128(%ebp)
#define  ihi          -136(%ebp)
#define  klor         -144(%ebp)
#define  khir         -152(%ebp)
#define  kloi         -160(%ebp)
#define  khii         -168(%ebp)
#define  U0ptr        -1192(%ebp)
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

         .global      butterfly4_and_transpose

         .data
         .align       32
t45:     .double      M_SQRT1_2
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
         mov          N2, %rbx
         mov          NMID, %rcx
         mov          NMID, %rdx
         mulx         %rbx, %rdx, %rdx
         shr          $2, %rax
         shr          $2, %rbx
         shl          $2, %rcx
         shl          $2, %rdx
         mov          %rax, NHIN1
         mov          %rbx, N2N1
         mov          %rcx, N1NMID
         mov          %rdx, N1N2NMID
         mov          N, %rax
         xor          %rdx, %rdx
         mov          N2, %rbx
         div          %rbx
         shr          %rax
         mov          %rax, TWHI
         xor          %rdx, %rdx
         xor          %rdx, ihi
ihiloop:
         xor          imid, imid
imidloop0:
         mov          N1NMID, %rax
         mulx         ihi, %rax, %rdi
         add          imid, %rdi
         mulx         N2, %rdi, %rdi
         shl          $2, %rdi
         xor          n0, n0
n0loop2:
         xor          n1, n1
n1loop2:
         mulx         N1NMID, n1, %rdx
         add          n0, %rdx
         mulx         N2, %rdx, %rdx
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
         vmovupd      (U0), er0
         vmovupd      32(U0), er1
         vmovupd      64(U0), er2
         vmovupd      96(U0), er3
         vaddpd       er2, er0, tr0
         vsubpd       er2, er0, tr2
         vaddpd       er3, er1, tr1
         vsubpd       er3, er1, tr3
         vaddpd       tr1, tr0, er0
         vsubpd       tr1, tr0, er2
         mulx         N1N2NMID, n0, %rax
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
         mov          NHIN1, %rax
         mulx         klor, %rax, %rdi
         add          ihi, %rdi
         mulx         N1NMID, %rdi, %rdi
         add          imid, %rdi
         mulx         N2, %rdi, %rdi
         add          khir, %rdi
         shl          $2, %rdi
         xor          n0, n0
n0loop4:
         xor          n1, n1
n1loop4:
         mulx         N1NMID, n1, %rdx
         add          n0, %rdx
         mulx         N2, %rdx, %rdx
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
n0loop5:
         mov          n0, %rdx
         shl          $7, %rdx
         add          U0, %rdx
         vmovupd      (U0), er0
         vmovupd      32(U0), er1
         vmovupd      64(U0), er2
         vmovupd      96(U0), er3
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
         mulx         N1N2NMID, n0, %rax
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
         mulx         TWHI, k2, %rax
         shl          %rax
         mulx         two, %rax, %rbx
         mulx         three, %rax, %rcx
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
         mov          NHIN1, %rax
         mulx         klor, %rax, %rdi
         mulx         kloi, %rax, %rsi
         add          ihi, %rdi
         add          ihi, %rsi
         mulx         N1NMID, %rdi, %rdi
         mulx         N1NMID, %rsi, %rsi
         add          imid, %rdi
         add          imid, %rsi
         mulx         N2, %rdi, %rdi
         mulx         N2, %rsi, %rsi
         add          khir, %rdi
         add          khii, %rsi
         shl          $2, %rdi
         shl          $2, %rsi
         xor          n0, n0
n0loop0:
         xor          n1, n1
n1loop0:
         mulx         N1NMID, n1, %rdx
         add          n0, %rdx
         mulx         N2, %rdx, %rdx
         mov          %rdx, %rcx
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
         vmovupd      (U0), er0
         vmovupd      32(U0), ei0
         vmovupd      64(U0), er1
         vmovupd      96(U0), ei1
         vmovupd      128(U0), er2
         vmovupd      160(U0), ei2
         vmovupd      192(U0), er3
         vmovupd      224(U0), ei3
         vbroadcastsd cos1, tr1
         vbroadcastsd cos2, tr2
         vbroadcastsd cos3, tr3
         vbroadcastsd sin1, ti1
         vbroadcastsd sin2, ti2
         vbroadcastsd sin3, ti3
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
         mulx         N1N2NMID, n0, %rdx
         mov          %rdx, %rcx
         add          %rdi, %rcx
         add          %rsi, %rdx
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
         cmp          %rax, NHIN1
         jne          ihiloop
         pop          %r15
         pop          %r14
         pop          %r13
         pop          %r12
         pop          %rbx
         leave
         ret

