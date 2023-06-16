#define  M_SQRT1_2    0.70710678118654752440
#define  N1           $4
#define  ilo          %r8
#define  k2           %r9
#define  depth        %r13
#define  X            %r11
#define  N            %r12
#define  NLO          %r10
#define  C            %r14
#define  S            %r15
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
#define  cos1         -8(%rbp)
#define  cos2         -16(%rbp)
#define  cos3         -24(%rbp)
#define  sin1         -32(%rbp)
#define  sin2         -40(%rbp)
#define  sin3         -48(%rbp)
#define  N2           -56(%rbp)
#define  STACK_SIZE   $56

         .global      fft_recursive

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

fft_recursive:
         enter        STACK_SIZE, $0
         push         %rbx
         push         %r12
         push         %r13
         push         %r14
         push         %r15
         mov          %rdi, X
         mov          %rsi, N
         mov          %rdx, C
         mov          %rcx, S
         mov          N, NLO
         shr          $2, NLO
         mov          $1, %rax
         mov          %rax, N2
         xor          depth, depth
         xor          k2, k2
         call         recurse
         pop          %r15
         pop          %r14
         pop          %r13
         pop          %r12
         pop          %rbx
         leave
         ret
recurse: cmp          N1, NLO
         jle          done
         push         k2
         mov          k2, %rdx
         mov          depth, %rcx
         shr          $2, %rcx
         mov          $1, %rax
         mov          $1, %rbx
         shl          %rcx, %rax
         inc          %rcx
         shl          %rcx, %rax
         or           %rax, k2
         or           %rbx, k2
         push         X
         push         k2
         not          %rax
         and          %rax, k2
         push         k2
         not          %rax
         not          %rbx
         or           %rax, k2
         and          %rbx, k2
         push         k2
         mov          %rdx, k2
         mov          N2, %rax
         shl          $2, %rax
         mov          N2, %rax
         inc          depth
         shr          $2, NLO
         call         recurse
         pop          k2
         lea          (X, NLO, 2), X
         call         recurse
         pop          k2
         lea          (X, NLO, 2), X
         call         recurse
         pop          k2
         lea          (X, NLO, 2), X
         call         recurse
         shl          $2, NLO
         mov          N2, %rax
         shr          $2, %rax
         mov          N2, %rax
         dec          depth
         pop          k2
         pop          X
         mov          N2, %rax
         shr          %rax
         cmp          k2, %rax
         jl           done
         je           N2o2
         cmp          $0, k2
         jne          k2rest
         xor          ilo, ilo
k0loop:  lea          (X, ilo, 8), %rax
         lea          (%rax, NLO, 2), %rbx
         lea          (%rbx, NLO, 2), %rcx
         lea          (%rcx, NLO, 2), %rdx
         vmovapd      (%rax), er0
         vmovapd      (%rbx), er2
         vmovapd      (%rcx), er1
         vmovapd      (%rdx), er3
         vaddpd       er2, er0, tr0
         vsubpd       er2, er0, tr2
         vaddpd       er3, er1, tr1
         vsubpd       er1, er3, tr3
         vaddpd       tr1, tr0, er0
         vsubpd       tr1, tr0, er2
         vmovapd      er0, (%rax)
         vmovapd      er2, (%rbx)
         vmovapd      tr2, (%rcx)
         vmovapd      tr3, (%rdx)
N2o2:
k2rest:  mov          k2, %rax
         imul         NLO, %rax
         imul         $2, %rax, %rbx
         imul         $3, %rax, %rcx
         movsd        (C, %rax, 8), %xmm0
         movsd        (C, %rbx, 8), %xmm1
         movsd        (C, %rcx, 8), %xmm2
         movsd        (S, %rax, 8), %xmm3
         movsd        (S, %rbx, 8), %xmm4
         movsd        (S, %rcx, 8), %xmm5
         movsd        %xmm0, cos1
         movsd        %xmm1, cos2
         movsd        %xmm2, cos3
         movsd        %xmm3, sin1
         movsd        %xmm4, sin2
         movsd        %xmm5, sin3
         push         %r15
         push         %r14
         push         %r13
         push         %r12
         xor          ilo, ilo
k2loop:  vbroadcastsd cos1, tr1
         vbroadcastsd cos2, tr2
         vbroadcastsd cos3, tr3
         vbroadcastsd sin1, ti1
         vbroadcastsd sin2, ti2
         vbroadcastsd sin3, ti3
         lea          (X, ilo, 8), %rax
         lea          (%rax, NLO, 2), %rbx
         lea          (%rbx, NLO, 2), %rcx
         lea          (%rcx, NLO, 2), %rdx
         lea          (%rdx, NLO, 2), %r12
         lea          (%r12, NLO, 2), %r13
         lea          (%r13, NLO, 2), %r14
         lea          (%r14, NLO, 2), %r15
         vmovapd      (%rax), er0
         vmovapd      (%rbx), er1
         vmovapd      (%rcx), er2
         vmovapd      (%rdx), er3
         vmovapd      (%r12), ei0
         vmovapd      (%r13), ei1
         vmovapd      (%r14), ei2
         vmovapd      (%r15), ei3
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
         vmovapd      ei3, (%rbx)
         vmovapd      er2, (%rcx)
         vmovapd      ei1, (%rdx)
         vmovapd      er1, (%r12)
         vmovapd      ei2, (%r13)
         vmovapd      er3, (%r14)
         vmovapd      ei0, (%r15)
         add          $4, ilo
         cmp          ilo, NLO
         jg           k2loop
done:    pop          %r12
         pop          %r13
         pop          %r14
         pop          %r15
         ret



