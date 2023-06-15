#define  M_SQRT1_2    0.70710678118654752440
#define  N1           $4
#define  X            %r8
#define  n0           %r9
#define  U0           %r10
#define  C            %r11
#define  S            %r12
#define  k2hi         %r13
#define  k2lo         %r14
#define  k2           %r15
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
#define  logNHINMIDN2 -56(%rbp)
#define  NHIN1N2NMID  -64(%rbp)


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
#define  N2           -176(%rbp)
#define  logN1N2NMID  -184(%rbp)
#define  logN2        -192(%rbp)
#define  imid         -200(%rbp)
#define  N1N2NMIDN1   -208(%rbp)
#define  dindex       -216(%rbp)
#define  U0ptr        -984(%rbp)
#define  STACK_SIZE   $1008

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
         shl          $2, %rax
         mov          %rax, N1N2NMIDN1
         shr          $2, %rax
         bsf          %rax, %rbx
         mov          %rbx, logN1N2NMID
         shr          $2, %rax
         mov          NHI, %rbx
         imul         %rbx, %rax
         bsf          %rax, %rbx
         mov          %rbx, logNHINMIDN2
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
         bsf          N2, %rdx
         mov          %rdx, logN2
         xor          %rax, %rax
         mov          %rax, ihi
ihi_loop0:
         xor          %rax, %rax
         mov          %rax, imid
imid_loop0:
         mov          imid, %rax
         add          ihi, %rax
         mov          %rax, dindex
         xor          n0, n0
k0_load:
         mov          n0, %rax
         mov          logN2, %rcx
         shl          %rcx, %rax
         lea          (X, %rax, 8), %rax
         mov          N1N2NMID, %rsi
         mov          dindex, %rdx
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
         mov          logN1N2NMID, %rcx
         mov          n0, %rax
         shl          %rcx, %rax
         add          dindex, %rax
         lea          (X, %rax, 8), %rax
         mov          N2, %rdx
         lea          (%rax, %rdx, 8), %rbx
         lea          (%rbx, %rdx, 8), %rcx
         lea          (%rcx, %rdx, 8), %rdx
         vmovapd      er0, (%rax)
         vmovapd      tr2, (%rbx)
         vmovapd      er2, (%rcx)
         vmovapd      tr3, (%rdx)
         dec          n0
         test         n0, n0
         jge          k0_compute
         mov          imid, %rax
         add          N1N2, %rax
         mov          %rax, imid
         cmp          %rax, N1N2NMID
         jg           imid_loop0
         mov          ihi, %rax
         add          N1N2NMIDN1, %rax
         mov          %rax, ihi
         cmp          %rax, NHIN1N2NMID
         jg           ihi_loop0
         xor          %rax, %rax
         mov          %rax, ihi
ihi_loop1:
         xor          %rax, %rax
         mov          %rax, imid
imid_loop1:
         mov          imid, %rax
         add          ihi, %rax
         mov          %rax, dindex
         xor          n0, n0
         mov          dN2o2, %rdi
         add          dindex, %rdi
         xor          n0, n0
kN2o2_load:
         mov          n0, %rax
         mov          logN2, %rcx
         shl          %rcx, %rax
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
         mov          logN1N2NMID, %rcx
         mov          n0, %rax
         shl          %rcx, %rax
         add          %rdi, %rax
         lea          (X, %rax, 8), %rax
         mov          N2, %rdx
         lea          (%rax, %rdx, 8), %rbx
         lea          (%rbx, %rdx, 8), %rcx
         lea          (%rcx, %rdx, 8), %rdx
         vmovapd      er0, (%rax)
         vmovapd      er1, (%rbx)
         vmovapd      er2, (%rcx)
         vmovapd      er3, (%rdx)
         dec          n0
         test         n0, n0
         jge          kN2o2_compute
         mov          imid, %rax
         add          N1N2, %rax
         mov          %rax, imid
         cmp          %rax, N1N2NMID
         jg           imid_loop1
         mov          ihi, %rax
         add          N1N2NMIDN1, %rax
         mov          %rax, ihi
         cmp          %rax, NHIN1N2NMID
         jg           ihi_loop1
         xor          k2lo, k2lo
k2lo_loop_b:
         cmp          N1, k2lo
         jge          k2lo_loop_e
         xor          %rax, %rax
         mov          %rax, ihi
ihi_loop2:
         xor          %rax, %rax
         mov          %rax, imid
imid_loop2:
         mov          imid, %rax
         add          ihi, %rax
         mov          %rax, dindex
         xor          k2hi, k2hi
k2hi_loop_b:
         mov          N2, %rax
         shr          $2, %rax
         cmp          k2hi, %rax
         je           k2hi_loop_e
         mov          k2hi, k2
         shl          $2, k2
         or           k2lo, k2
         cmp          $0, k2
         je           k2hi_loop_c
         cmp          N2o2, k2
         jge          k2hi_loop_c
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
         mov          logNHINMIDN2, %rcx
         mov          k2lo, %rdi
         shl          %rcx, %rdi
         add          k2hi, %rdi
         shl          $2, %rdi
         mov          N2, %rsi
         sub          k2, %rsi
         mov          %rsi, %rax
         shr          $2, %rax
         and          $3, %rsi
         shl          %rcx, %rsi
         add          %rax, %rsi
         shl          $2, %rsi
         sub          %rdi, %rsi
         add          dindex, %rdi




         xor          n0, n0
load_loop:
         mov          n0, %rax
         imul         N2, %rax
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
         mov          logN1N2NMID, %rcx
         mov          n0, %rax
         shl          %rcx, %rax
         add          %rdi, %rax
         lea          (X, %rax, 8), %rax
         mov          N2, %rdx
         lea          (%rax, %rdx, 8), %rbx
         lea          (%rbx, %rdx, 8), %rcx
         lea          (%rcx, %rdx, 8), %rdx
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












k2hi_loop_c:
         inc          k2hi
         jmp          k2hi_loop_b
k2hi_loop_e:
         mov          imid, %rax
         add          N1N2, %rax
         mov          %rax, imid
         cmp          %rax, N1N2NMID
         jg           imid_loop2
         mov          ihi, %rax
         add          N1N2NMIDN1, %rax
         mov          %rax, ihi
         cmp          %rax, NHIN1N2NMID
         jg           ihi_loop2
         inc          k2lo
         jmp          k2lo_loop_b
k2lo_loop_e:
         pop          %r15
         pop          %r14
         pop          %r13
         pop          %r12
         pop          %rbx
         leave
         ret

