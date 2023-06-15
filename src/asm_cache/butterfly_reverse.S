#define  M_SQRT1_2    0.70710678118654752440
#define  N2           %r8
#define  k2rev        %r9
#define  k2           %r10
#define  nlo          %r11
#define  NLO          %r12
#define  C            %r13
#define  S            %r14
#define  X            %r15
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
#define  TW1          -8(%rbp)
#define  TW2          -16(%rbp)
#define  TW3          -24(%rbp)
#define  cos1         -32(%rbp)
#define  cos2         -40(%rbp)
#define  cos3         -48(%rbp)
#define  sin1         -56(%rbp)
#define  sin2         -64(%rbp)
#define  sin3         -72(%rbp)
#define  N2o2         -80(%rbp)
#define  No2          -88(%rbp)
#define  N1NLO        -96(%rbp)
#define  N2m1         -104(%rbp)
#define  N2NLOo2      -112(%rbp)
#define  STACK_SIZE   $112

         .global      butterfly_reverse

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

butterfly_reverse:
         enter        STACK_SIZE, $0
         push         %r15
         push         %r14
         push         %r13
         push         %r12
         push         %rbx
         mov          %r8, NLO
         mov          %rcx, N2
         mov          %rdx, S
         mov          %rsi, C
         mov          %rdi, X
         mov          N2, %rbx
         mov          N2, %rax
         dec          %rbx
         mov          %rbx, N2m1
         shl          %rax
         mov          %rax, No2
         shr          $2, %rax
         mov          %rax, N2o2
         mov          NLO, %rax
         mov          NLO, %rbx
         mov          NLO, %rcx
         mov          NLO, %rdx
         shl          $2, %rdx
         mov          %rdx, N1NLO
         mov          N2, %rdx
         imul         NLO, %rdx
         shr          %rdx
         mov          %rdx, N2NLOo2
         add          %rax, %rbx
         add          %rbx, %rax
         mov          %rax, TW1
         mov          %rbx, TW2
         mov          %rcx, TW3
         xor          nlo, nlo
nlo_loop0:
         lea          (X, nlo, 8), %rax
         lea          (%rax, NLO, 8), %rbx
         lea          (%rbx, NLO, 8), %rcx
         lea          (%rcx, NLO, 8), %rdx
         vmovapd      (%rax), er0
         vmovapd      (%rbx), er1
         vmovapd      (%rcx), er2
         vmovapd      (%rdx), er3
         vaddpd       er2, er0, tr0
         vsubpd       er2, er0, tr2
         vaddpd       er3, er1, tr1
         vsubpd       er1, er3, tr3
         vaddpd       tr1, tr0, er0
         vsubpd       tr1, tr0, er2
         vmovapd      er0, (%rax)
         vmovapd      tr2, (%rbx)
         vmovapd      er2, (%rcx)
         vmovapd      tr3, (%rdx)
         add          $4, nlo
         cmp          nlo, NLO
         jg           nlo_loop0
         xor          nlo, nlo
nlo_loop1:
         mov          N2NLOo2, %rdi
         add          nlo, %rdi
         lea          (X, %rdi, 8), %rax
         lea          (%rax, NLO, 8), %rbx
         lea          (%rbx, NLO, 8), %rcx
         lea          (%rcx, NLO, 8), %rdx
         vmovapd      (%rax), er0
         vmovapd      (%rbx), er1
         vmovapd      (%rcx), er2
         vmovapd      (%rdx), er3
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
         vmovapd      er0, (%rax)
         vmovapd      er1, (%rbx)
         vmovapd      er2, (%rcx)
         vmovapd      er3, (%rdx)
         add          $4, nlo
         cmp          nlo, NLO
         jg           nlo_loop1
         mov          $1, k2
         mov          No2, k2rev
k2_begin:
         mov          k2rev, %rdx
         mulx         TW1, %rax, %rcx
         mulx         TW2, %rbx, %rsi
         mulx         TW3, %rcx, %rdi
         movq         (C, %rax, 8), %rdx
         movq         (S, %rax, 8), %rax
         movq         (C, %rbx, 8), %rdi
         movq         (S, %rbx, 8), %rbx
         movq         (C, %rcx, 8), %rsi
         movq         (S, %rcx, 8), %rcx
         movq         %rdx, cos1
         movq         %rax, sin1
         movq         %rdi, cos2
         movq         %rbx, sin2
         movq         %rsi, cos3
         movq         %rcx, sin3
         xor          nlo, nlo
nlo_loop2:
         vbroadcastsd cos1, tr1
         vbroadcastsd cos2, tr2
         vbroadcastsd cos3, tr3
         vbroadcastsd sin1, ti1
         vbroadcastsd sin2, ti2
         vbroadcastsd sin3, ti3
         mov          N2, %rsi
         sub          k2, %rsi
         mov          k2, %rdi
         sub          k2, %rsi
         imul         N1NLO, %rdi
         imul         N1NLO, %rsi
         add          nlo, %rdi
         lea          (X, %rdi, 8), %rax
         lea          (%rax, NLO, 8), %rbx
         lea          (%rbx, NLO, 8), %rcx
         lea          (%rcx, NLO, 8), %rdx
         vmovapd      (%rax), er0
         vmovapd      (%rbx), er1
         vmovapd      (%rcx), er2
         vmovapd      (%rdx), er3
         vmovapd      (%rax, %rsi, 8), ei0
         vmovapd      (%rbx, %rsi, 8), ei1
         vmovapd      (%rcx, %rsi, 8), ei2
         vmovapd      (%rdx, %rsi, 8), ei3
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
         vmovapd      er0, (%rax)
         vmovapd      er1, (%rbx)
         vmovapd      ei2, (%rcx)
         vmovapd      ei3, (%rdx)
         vmovapd      er3, (%rax, %rsi, 8)
         vmovapd      er2, (%rbx, %rsi, 8)
         vmovapd      ei1, (%rcx, %rsi, 8)
         vmovapd      ei0, (%rdx, %rsi, 8)
         add          $4, nlo
         cmp          nlo, NLO
         jg           nlo_loop2
         cmp          k2, N2m1
         je           k2_end
L80:     mov          No2, %rdx
L90:     cmp          %rdx, k2rev
         jl           L100
         sub          %rdx, k2rev
         shr          %rdx
         jmp          L90
L100:    add          %rdx, k2rev
L110:    inc          k2
         jmp          k2_begin
k2_end:
         pop          %rbx
         pop          %r12
         pop          %r13
         pop          %r14
         pop          %r15
         leave
         ret


