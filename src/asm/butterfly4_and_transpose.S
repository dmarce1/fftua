#define  M_SQRT1_2    0.70710678118654752440
#define  N1           $4
#define  X            %r8
#define  C            %r9
#define  S            %r10
#define  N2           %r11
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
#define  N            -8(%rbp)
#define  N2o2         -16(%rbp)
#define  NHI          -24(%rbp)
#define  N1N2         -32(%rbp)
#define  NMID         -40(%rbp)
#define  N1N2NMID     -48(%rbp)
#define  NHINMIDN2    -56(%rbp)
#define  NHIN1N2NMID  -64(%rbp)
#define  mask1        -72(%rbp)
#define  mask2        -80(%rbp)
#define  cos1         -88(%rbp)
#define  cos2         -96(%rbp)
#define  cos3         -104(%rbp)
#define  sin1         -112(%rbp)
#define  sin2         -120(%rbp)
#define  sin3         -128(%rbp)
#define  TW1          -136(%rbp)
#define  TW2          -144(%rbp)
#define  TW3          -152(%rbp)
#define  dN2o2        -160(%rbp)
#define  ihi          -168(%rbp)
#define  U0ptr        -936(%rbp)
#define  STACK_SIZE   $960

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
mask_lo: .quad        512
mask_hi: .quad        7426


         .text

butterfly4_and_transpose:
         enter        STACK_SIZE, $0
         push         %rbx
         push         %r12
         push         %r13
         push         %r14
         push         %r15
         mov          %r9, N2
         mov          %r8, NMID
         mov          %rcx, NHI
         mov          %rdx, S
         mov          %rsi, C
         mov          %rdi, X
         lea          U0ptr, U0
         and          $0xffffffffffffffe0, U0
         mov          N2, %rax
         shl          $2, %rax
         mov          %rax, N1N2
         mov          NMID, %rbx
         imul         %rbx, %rax
         mov          %rax, N1N2NMID
         shr          $2, %rax
         mov          NHI, %rbx
         imul         %rbx, %rax
         mov          %rax, NHINMIDN2
         shl          $2, %rax
         mov          %rax, NHIN1N2NMID
         mov          N2, %rbx
         shr          %rbx
         mov          %rbx, N2o2
         mov          NMID, %rax
         shl          $2, %rax
         imul         NHI, %rax
         shl          $3, %rax
         mov          %rax, %rbx
         mov          %rax, TW1
         add          %rax, %rbx
         mov          %rbx, TW2
         add          %rax, %rbx
         mov          %rbx, TW3
         mov          NMID, %rdx
         bsf          %rdx, %rcx
         mov          N2, %rdx
         bsf          %rdx, %rbx
         add          %rbx, %rcx
         add          $2, %rcx
         mov          $3, %rax
         shl          %cl, %rax
         mov          %rax, mask1
         not          %rax
         mov          %rax, mask2
         mov          NHI, %rdi
         mov          N2o2, %rax
         and          $3, %rax
         imul         %rax, %rdi
         imul         NMID, %rdi
         imul         N2, %rdi
         mov          N2o2, %rax
         shr          $2, %rax
         add          %rax, %rdi
         shl          $2, %rdi
         mov          %rdi, dN2o2
         xor          %rax, %rax
         mov          %rax, ihi
ihi_loop:
         xor          n0, n0
k0_load:
         mov          n0, %rax
         imul         N2
         lea          (X, %rax, 8), %rax
         mov          N1N2NMID, %rsi
         mov          ihi, %rdx
         lea          (%rax, %rdx, 8), %rax
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
         add          ihi, %rax
         lea          (X, %rax, 8), %rax
         lea          (%rax, N2, 8), %rbx
         lea          (%rbx, N2, 8), %rcx
         lea          (%rcx, N2, 8), %rdx
         vmovapd      er0, (%rax)
         vmovapd      tr2, (%rbx)
         vmovapd      er2, (%rcx)
         vmovapd      tr3, (%rdx)
         dec          n0
         test         n0, n0
         jge          k0_compute
         mov          dN2o2, %rdi
         add          ihi, %rdi
         xor          n0, n0
kN2o2_load:
         mov          n0, %rax
         imul         N2, %rax
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
         je           kN2o2_enter
         mov          n0, %rax
         shl          $4, %rax
         lea          (U0, %rax, 8), %rax
         vmovapd      er0, 0(%rax)
         vmovapd      er1, 32(%rax)
         vmovapd      er2, 64(%rax)
         vmovapd      er3, 96(%rax)
         inc          n0
         jmp          kN2o2_load
kN2o2_compute:
         mov          n0, %rdx
         shl          $7, %rdx
         add          U0, %rdx
         vmovapd      (%rdx), er0
         vmovapd      32(%rdx), er1
         vmovapd      64(%rdx), er2
         vmovapd      96(%rdx), er3
kN2o2_enter:
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
         lea          (%rax, N2, 8), %rbx
         lea          (%rbx, N2, 8), %rcx
         lea          (%rcx, N2, 8), %rdx
         vmovapd      er0, (%rax)
         vmovapd      er1, (%rbx)
         vmovapd      er2, (%rcx)
         vmovapd      er3, (%rdx)
         dec          n0
         test         n0, n0
         jge          kN2o2_compute
         xor          k2, k2
k2_loop_b:
         inc          k2
         cmp          k2, N2o2
         je           k2_loop_e
         mov          k2, %rdx
         mulx         TW1, %rax, %rdi
         mulx         TW2, %rbx, %rdi
         mulx         TW3, %rcx, %rdi
         mov          (C, %rax), %rdi
         mov          (S, %rax), %rax
         mov          (C, %rbx), %rsi
         mov          (S, %rbx), %rbx
         mov          (C, %rcx), %rdx
         mov          (S, %rcx), %rcx
         mov          %rdi, cos1
         mov          %rax, sin1
         mov          %rsi, cos2
         mov          %rbx, sin2
         mov          %rdx, cos3
         mov          %rcx, sin3
         mov          mask_hi, %rax
         mov          mask_lo, %rsi
         mov          N2, %rbx
         sub          k2, %rbx
         bextr        %rax, k2, %rcx
         bextr        %rsi, k2, %rdi
         bextr        %rax, %rbx, %rax
         bextr        %rsi, %rbx, %rsi
         imul         NHINMIDN2, %rdi
         imul         NHINMIDN2, %rsi
         add          %rcx, %rdi
         add          %rax, %rsi
         shl          $2, %rsi
         shl          $2, %rdi
         sub          %rdi, %rsi
         add          ihi, %rdi
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
         vmovapd      (%rbx), %ymm2
         vmovapd      (%rcx), %ymm4
         vmovapd      (%rdx), %ymm6
         vmovapd      (%rax, %rsi, 8), %ymm1
         vmovapd      (%rbx, %rsi, 8), %ymm3
         vmovapd      (%rcx, %rsi, 8), %ymm5
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
         test         n0, n0
         jge          compute_loop
         jmp          k2_loop_b
k2_loop_e:
         mov          ihi, %rax
         or           mask1, %rax
         add          N1N2, %rax
         and          mask2, %rax
         mov          %rax, ihi
         cmp          %rax, NHIN1N2NMID
         jg           ihi_loop
         pop          %r15
         pop          %r14
         pop          %r13
         pop          %r12
         pop          %rbx
         leave
         ret

