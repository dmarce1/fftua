#define  STACK_SIZE   $1080
#define  M_SQRT1_2    0.70710678118654752440
#define  N1           $4
#define  X            %r8
#define  W            %r9
#define  N2           %r10
#define  imid         %r11
#define  k2           %r12
#define  n0           %r13
#define  n1           %r14
#define  U0           %r15
#define  er0          %ymm0
#define  ei0          %ymm1
#define  er1          %ymm2
#define  ei1          %ymm3
#define  er2          %ymm4
#define  ei2          %ymm5
#define  er3          %ymm6
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
#define  N1N2         -16(%rbp)
#define  NMID         -24(%rbp)
#define  NHIoN1       -32(%rbp)
#define  N1NMID       -40(%rbp)
#define  N1N2NMID     -48(%rbp)
#define  N            -56(%rbp)
#define  ihi          -64(%rbp)
#define  N2o2         -72(%rbp)
#define  cos1         -80(%rbp)
#define  cos2         -112(%rbp)
#define  cos3         -144(%rbp)
#define  sin1         -176(%rbp)
#define  sin2         -208(%rbp)
#define  sin3         -240(%rbp)
#define  TW1          -248(%rbp)
#define  TW2          -256(%rbp)
#define  TW3          -264(%rbp)
#define  NHINMIDN2    -272(%rbp)
#define  N2o2hi       -280(%rbp)
#define  N2o2lo       -288(%rbp)
#define  U0ptr        -1056(%rbp)


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
         and          $0xffffffffffffffe0, U0
         mov          NHI, %rax
         shr          $2, %rax
         mov          %rax, NHIoN1
         mov          NMID, %rbx
         shl          $2, %rbx
         mov          %rbx, N1NMID
         mov          N2, %rax
         imul         %rbx, %rax
         mov          %rax, N1N2NMID
         mov          NHI, %rbx
         imul         %rbx, %rax
         shl          $2, %rax
         mov          %rax, N
         mov          N, %rax
         mov          N2, %rbx
         shl          $2, %rbx
         xor          %rdx, %rdx
         div          %rbx
         shl          $4, %rax
         mov          %rax, TW1
         mov          %rax, %rbx
         add          %rax, %rbx
         mov          %rbx, TW2
         add          %rax, %rbx
         mov          %rbx, TW3
         mov          N2, %rax
         shl          $2, %rax
         mov          %rax, N1N2
         shr          $3, %rax
         mov          %rax, N2o2
         xor          %rdx, %rdx
         mov          NHI, %rax
         mov          NMID, %rbx
         imul         %rbx, %rax
         imul         N2, %rax
         mov          %rax, NHINMIDN2
         mov          N2o2, %rax
         and          $3, %rax
         mov          %rax, N2o2lo
         mov          N2o2, %rax
         shr          $2, %rax
         mov          %rax, N2o2hi
         mov          %rdx, ihi
ihi_loop:
         xor          imid, imid
imid_k0:
         mov          N1NMID, %rax
         imul         ihi, %rax
         add          imid, %rax
         imul         N2, %rax
         shl          $2, %rax
         mov          %rax, %rdi
         xor          n0, n0
k0_load:
         mov          n0, %rax
         imul         N2
         lea          (X, %rax, 8), %rax
         mov          N1N2NMID, %rsi
         lea          (%rax, %rdi, 8), %rax
         lea          (%rax, %rsi, 8), %rbx
         lea          (%rbx, %rsi, 8), %rcx
         lea          (%rcx, %rsi, 8), %rdx
         vmovapd      (%rax), er0
         vmovapd      (%rbx), er1
         vmovapd      (%rcx), er2
         vmovapd      (%rdx), er3
         cmp          $3, n0
         je           k0_compute_enter
         mov          n0, %rax
         shl          $4, %rax
         lea          (U0, %rax, 8), %rax
         vmovapd      er0, 0(%rax)
         vmovapd      er1, 32(%rax)
         vmovapd      er2, 64(%rax)
         vmovapd      er3, 96(%rax)
         inc          n0
         jmp          k0_load
k0_compute:
         mov          n0, %rdx
         shl          $7, %rdx
         add          U0, %rdx
         vmovapd      (%rdx), er0
         vmovapd      32(%rdx), er1
         vmovapd      64(%rdx), er2
         vmovapd      96(%rdx), er3
k0_compute_enter:
         vaddpd       er2, er0, tr0
         vsubpd       er2, er0, tr2
         vaddpd       er3, er1, tr1
         vsubpd       er1, er3, tr3
         vaddpd       tr1, tr0, er0
         vsubpd       tr1, tr0, er2
         mov          N1N2NMID, %rax
         imul         n0, %rax
         add          %rdi, %rax
         lea          (X, %rax, 8), %rax
         lea          (%rax, N2, 8), %rbx
         lea          (%rbx, N2, 8), %rcx
         lea          (%rcx, N2, 8), %rdx
         vmovapd      er0, (%rax)
         vmovapd      tr2, (%rbx)
         vmovapd      er2, (%rcx)
         vmovapd      tr3, (%rdx)
         dec          n0
         cmp          $0, n0
         jge          k0_compute
         inc          imid
         cmp          imid, NMID
         jne          imid_k0
         xor          imid, imid
imid_kN2o2:
         mov          NHIoN1, %rdi
         imul         N2o2lo, %rdi
         add          ihi, %rdi
         imul         N1NMID, %rdi
         add          imid, %rdi
         imul         N2, %rdi
         add          N2o2hi, %rdi
         shl          $2, %rdi
         xor          n0, n0
kN2o2_load:
         mov          n0, %rdx
         mulx         N2, %rcx, %rax
         lea          (X, %rcx, 8), %rax
         lea          (%rax, %rdi, 8), %rcx
         mov          N1N2NMID, %rax
         vmovapd      (%rcx), %ymm0
         lea          (%rcx, %rax, 8), %rcx
         vmovapd      (%rcx), %ymm1
         lea          (%rcx, %rax, 8), %rcx
         vmovapd      (%rcx), %ymm2
         lea          (%rcx, %rax, 8), %rcx
         vmovapd      (%rcx), %ymm3
         lea          (%rcx, %rax, 8), %rcx
         mov          n0, %rax
         shl          $4, %rax
         lea          (U0, %rax, 8), %rax
         vmovapd      %ymm0, 0(%rax)
         vmovapd      %ymm1, 32(%rax)
         vmovapd      %ymm2, 64(%rax)
         vmovapd      %ymm3, 96(%rax)
         inc          n0
         cmp          N1, n0
         jne          kN2o2_load
         xor          n0, n0
kN2o2_compute:
         mov          n0, %rdx
         shl          $7, %rdx
         add          U0, %rdx
         vmovapd      (%rdx), er0
         vmovapd      32(%rdx), er1
         vmovapd      64(%rdx), er2
         vmovapd      96(%rdx), er3
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
         mul          n0
         add          %rdi, %rax
         lea          (X, %rax, 8), %rax
         vmovapd      er0, (%rax)
         lea          (%rax, N2, 8), %rax
         vmovapd      er1, (%rax)
         lea          (%rax, N2, 8), %rax
         vmovapd      er2, (%rax)
         lea          (%rax, N2, 8), %rax
         vmovapd      er3, (%rax)
         inc          n0
         cmp          N1, n0
         jne          kN2o2_compute
         inc          imid
         cmp          imid, NMID
         jne          imid_kN2o2
         xor          imid, imid
imid_loop:
         mov          ihi, %rax
         imul         N1NMID, %rax
         add          imid, %rax
         imul         N1N2, %rax
         push         %rax
         xor          k2, k2
k2_loop_b:
         inc          k2
         cmp          k2, N2o2
         je           k2_loop_e
         mov          k2, %rdx
         mulx         TW1, %rax, %rdi
         mulx         TW2, %rbx, %rdi
         mulx         TW3, %rcx, %rdi
         mov          0(W, %rax), %rdi
         mov          8(W, %rax), %rax
         mov          0(W, %rbx), %rsi
         mov          8(W, %rbx), %rbx
         mov          0(W, %rcx), %rdx
         mov          8(W, %rcx), %rcx
         mov          %rdi, cos1
         mov          %rax, sin1
         mov          %rsi, cos2
         mov          %rbx, sin2
         mov          %rdx, cos3
         mov          %rcx, sin3
         mov          N2, %rax
         mov          k2, %rcx
         sub          k2, %rax
         mov          k2, %rdx
         mov          %rax, %rbx
         and          $3, %rdx
         and          $3, %rbx
         shr          $2, %rcx
         shr          $2, %rax
         imul         NHINMIDN2, %rdx
         imul         NHINMIDN2, %rbx
         add          %rcx, %rdx
         add          %rbx, %rax
         shl          $2, %rdx
         shl          $2, %rax
         mov          %rdx, %rdi
         mov          %rax, %rsi
         pop          %rdx
         add          %rdx, %rdi
         add          %rdx, %rsi
         push         %rdx
         sub          %rdi, %rsi
         xor          n0, n0
load_loop:
         mov          n0, %rdx
         mulx         N2, %rax, %rcx
         add          %rdi, %rax
         mov          N1N2NMID, %rdx
         lea          (X, %rax, 8), %rax
         lea          (%rax, %rdx, 8), %rbx
         lea          (%rbx, %rdx, 8), %rcx
         lea          (%rcx, %rdx, 8), %rdx
         vmovapd      (%rax), %ymm0
         vmovapd      (%rax, %rsi, 8), %ymm1
         vmovapd      (%rbx), %ymm2
         vmovapd      (%rbx, %rsi, 8), %ymm3
         vmovapd      (%rcx), %ymm4
         vmovapd      (%rcx, %rsi, 8), %ymm5
         vmovapd      (%rdx), %ymm6
         vmovapd      (%rdx, %rsi, 8), %ymm7
         cmp          $3, n0
         je           compute_enter
         mov          n0, %rax
         shl          $5, %rax
         lea          (U0, %rax, 8), %rax
         vmovapd      %ymm0, 0(%rax)
         vmovapd      %ymm1, 32(%rax)
         vmovapd      %ymm2, 64(%rax)
         vmovapd      %ymm3, 96(%rax)
         vmovapd      %ymm4, 128(%rax)
         vmovapd      %ymm5, 160(%rax)
         vmovapd      %ymm6, 192(%rax)
         vmovapd      %ymm7, 224(%rax)
         inc          n0
         jmp          load_loop
compute_loop:
         mov          n0, %rdx
         shl          $8, %rdx
         add          U0, %rdx
         vmovapd      (%rdx), er0
         vmovapd      32(%rdx), ei0
         vmovapd      64(%rdx), er1
         vmovapd      96(%rdx), ei1
         vmovapd      128(%rdx), er2
         vmovapd      160(%rdx), ei2
         vmovapd      192(%rdx), er3
         vmovapd      224(%rdx), ei3
compute_enter:
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
         mov          n0, %rax
         imul         N1N2NMID, %rax
         add          %rdi, %rax
         lea          (X, %rax, 8), %rax
         lea          (%rax, N2, 8), %rbx
         lea          (%rbx, N2, 8), %rcx
         lea          (%rcx, N2, 8), %rdx
         vmovapd      er0, (%rax)
         vmovapd      er1, (%rbx)
         vmovapd      ei2, (%rcx)
         vmovapd      ei3, (%rdx)
         vmovapd      er3, (%rax, %rsi, 8)
         vmovapd      er2, (%rbx, %rsi, 8)
         vmovapd      ei1, (%rcx, %rsi, 8)
         vmovapd      ei0, (%rdx, %rsi, 8)
         dec          n0
         cmp          $0, n0
         jge          compute_loop
         jmp          k2_loop_b
k2_loop_e:
         pop          %rax
         inc          imid
         cmp          imid, NMID
         jne          imid_loop
         mov          ihi, %rax
         inc          %rax
         mov          %rax, ihi
         cmp          %rax, NHIoN1
         jne          ihi_loop
         pop          %r15
         pop          %r14
         pop          %r13
         pop          %r12
         pop          %rbx
         leave
         ret

