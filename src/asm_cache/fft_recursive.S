#define  M_SQRT1_2    0.70710678118654752440
#define  N1           $4
#define  X0           %rsi
#define  ilo          %r8
#define  k2           %r9
#define  NLO          %r10
#define  N2           %r11
#define  N            %r12
#define  X            %r13
#define  C            %r14
#define  S            %r15
#define  er0          %ymm0
#define  er1          %ymm1
#define  er2          %ymm2
#define  er3          %ymm3
#define  ei0          %ymm4
#define  ei1          %ymm5
#define  ei2          %ymm6
#define  ei3          %ymm7
#define  tr0          %ymm8
#define  ti0          %ymm9
#define  cos1         %ymm10
#define  cos2         %ymm11
#define  sin1         %ymm12
#define  sin2         %ymm13
#define  tr1          %ymm10
#define  tr2          %ymm11
#define  ti1          %ymm12
#define  ti2          %ymm13
#define  tr3          %ymm14
#define  ytwo         %ymm15
#define  ur0          %xmm0
#define  ui0          %xmm1
#define  ur1          %xmm2
#define  ui1          %xmm3
#define  ur2          %xmm4
#define  ui2          %xmm5
#define  ur3          %xmm6
#define  ui3          %xmm7
#define  sr0          %xmm8
#define  si0          %xmm9
#define  tcos1        %xmm10
#define  tcos2        %xmm11
#define  tsin1        %xmm12
#define  tsin2        %xmm13
#define  sr1          %xmm10
#define  sr2          %xmm11
#define  si1          %xmm12
#define  si2          %xmm13
#define  sr3          %xmm14
#define  xtwo         %xmm15
#define  prmt_cntrl   $27

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
two:     .double      2.0
         .double      2.0
         .double      2.0
         .double      2.0
pbytes:  .byte        8
         .byte        9
         .byte        10
         .byte        11
         .byte        12
         .byte        13
         .byte        14
         .byte        15
         .byte        0
         .byte        1
         .byte        2
         .byte        3
         .byte        4
         .byte        5
         .byte        6
         .byte        7


         .text

fft_recursive:
         push         %rbx
         push         %r12
         push         %r13
         push         %r14
         push         %r15
         mov          %rdi, X
         mov          %rsi, C
         mov          %rdx, S
         mov          %rcx, N
         mov          N, NLO
         shr          $2, NLO
         mov          $1, N2
         xor          k2, k2
         mov          X, X0
         call         recurse
         mov          X, %rdi
         mov          N, %rsi
         push         %r8
         push         %r9
         push         %r10
         push         %r11
         call         fft_scramble
         pop          %r11
         pop          %r10
         pop          %r9
         pop          %r8
         mov          N, N2
         shr          $2, N2
         lea          (X), %rax
         lea          (%rax, N2, 8), %rbx
         lea          (%rbx, N2, 8), %rcx
         lea          (%rcx, N2, 8), %rdx
         movsd        (%rax), ur0
         movsd        (%rbx), ur2
         movsd        (%rcx), ur1
         movsd        (%rdx), ur3
         vaddsd       ur2, ur0, sr0
         vsubsd       ur2, ur0, sr2
         vaddsd       ur3, ur1, sr1
         vsubsd       ur1, ur3, sr3
         vaddsd       sr1, sr0, ur0
         vsubsd       sr1, sr0, ur2
         movsd        ur0, (%rax)
         movsd        tcos2, (%rbx)
         movsd        ur2, (%rcx)
         movsd        sr3, (%rdx)
         cmp          $1, N2
         jle          done0
         movsd        (%rax, N), ur0
         movsd        (%rbx, N), ur2
         movsd        (%rcx, N), ur1
         movsd        (%rdx, N), ur3
         vaddsd       ur3, ur1, sr0
         vsubsd       ur3, ur1, sr2
         vmulsd       tw45, sr2, sr1
         vmulsd       tw45, sr0, sr3
         movsd        ur0, sr0
         movsd        ur2, tcos2
         vaddsd       sr1, sr0, ur0
         vaddsd       sr3, sr2, ur3
         vsubsd       sr1, sr0, ur1
         vsubsd       sr3, sr2, ur2
         vmulsd       none, ur3, ur3
         movsd        ur0, (%rax, N)
         movsd        ur1, (%rbx, N)
         movsd        ur2, (%rcx, N)
         movsd        ur3, (%rdx, N)
         cmp          $2, N2
         jle          done0
         movsd        8(C), tcos1
         movsd        8(S), tsin1
         vmulsd       tsin1, tsin1, tcos2
         vmulsd       tcos1, tsin1, tsin2
         vfmsub231sd  tcos1, tcos1, tcos2
         vfmadd231sd  tsin1, tcos1, tsin2
         mov          N2, %rdi
         sub          $2, %rdi
         lea          8(X), %rax
         lea          (%rax, N2, 8), %rbx
         lea          (%rbx, N2, 8), %rcx
         lea          (%rcx, N2, 8), %rdx
         movsd        (%rax), ur0
         movsd        (%rbx), ur2
         movsd        (%rcx), ur1
         movsd        (%rdx), ur3
         movsd        (%rax, %rdi, 8), ui0
         movsd        (%rbx, %rdi, 8), ui2
         movsd        (%rcx, %rdi, 8), ui1
         movsd        (%rdx, %rdi, 8), ui3
         movsd        ur0, sr0
         movsd        ui0, si0
         vmovapd      two, ytwo
         vfmadd231sd  tsin2, ui2, sr0
         vfnmadd231sd tsin2, ur2, si0
         vfnmadd132sd tcos2, sr0, ur2
         vfnmadd132sd tcos2, si0, ui2
         vfmsub132sd  xtwo, ur2, ur0
         vfmsub132sd  xtwo, ui2, ui0
         movsd        ur1, sr0
         movsd        ui1, si0
         vfmadd231sd  tsin2, ui3, sr0
         vfnmadd231sd tsin2, ur3, si0
         vfnmadd132sd tcos2, sr0, ur3
         vfnmadd132sd tcos2, si0, ui3
         vfmsub132sd  xtwo, ur3, ur1
         vfmsub132sd  xtwo, ui3, ui1
         movsd        ur0, sr0
         movsd        ui0, si0
         vfmadd231sd  tsin1, ui1, sr0
         vfnmadd231sd tsin1, ur1, si0
         vfnmadd132sd tcos1, sr0, ur1
         vfmsub132sd  tcos1, si0, ui1
         vfmsub132sd  xtwo, ur1, ur0
         vfmadd132sd  xtwo, ui1, ui0
         movsd        ur2, sr0
         movsd        ui2, si0
         vfmadd231sd  tcos1, ui3, sr0
         vfnmadd231sd tcos1, ur3, si0
         vfmadd132sd  tsin1, sr0, ur3
         vfmadd132sd  tsin1, si0, ui3
         vfmsub132sd  xtwo, ur3, ur2
         vfnmadd132sd xtwo, ui3, ui2
         movsd        ur0, (%rax)
         movsd        ur3, (%rbx)
         movsd        ui1, (%rcx)
         movsd        ui2, (%rdx)
         movsd        ur2, (%rax, %rdi, 8)
         movsd        ur1, (%rbx, %rdi, 8)
         movsd        ui3, (%rcx, %rdi, 8)
         movsd        ui0, (%rdx, %rdi, 8)
         cmp          $4, N2
         jle          done0
         movlpd       16(C), tcos1
         movlpd       16(S), tsin1
         vmovhpd      24(C), tcos1, tcos1
         vmovhpd      24(S), tsin1, tsin1
         vmulpd       sin1, sin1, cos2
         vmulpd       cos1, sin1, sin2
         vfmsub231pd  cos1, cos1, cos2
         vfmadd231pd  sin1, cos1, sin2
         mov          N2, %rdi
         sub          $5, %rdi
         lea          16(X), %rax
         lea          (%rax, N2, 8), %rbx
         lea          (%rbx, N2, 8), %rcx
         lea          (%rcx, N2, 8), %rdx
         vmovapd      (%rax), ur0
         vmovapd      (%rbx), ur2
         vmovapd      (%rcx), ur1
         vmovapd      (%rdx), ur3
         vmovupd      (%rax, %rdi, 8), ui0
         vmovupd      (%rbx, %rdi, 8), ui2
         vmovupd      (%rcx, %rdi, 8), ui1
         vmovupd      (%rdx, %rdi, 8), ui3
         pshufb       pbytes, ui0
         pshufb       pbytes, ui1
         pshufb       pbytes, ui2
         pshufb       pbytes, ui3
         vmovapd      ur0, sr0
         vmovapd      ui0, si0
         vfmadd231pd  tsin2, ui2, sr0
         vfnmadd231pd tsin2, ur2, si0
         vfnmadd132pd tcos2, sr0, ur2
         vfnmadd132pd tcos2, si0, ui2
         vfmsub132pd  xtwo, ur2, ur0
         vfmsub132pd  xtwo, ui2, ui0
         vmovapd      ur1, sr0
         vmovapd      ui1, si0
         vfmadd231pd  tsin2, ui3, sr0
         vfnmadd231pd tsin2, ur3, si0
         vfnmadd132pd tcos2, sr0, ur3
         vfnmadd132pd tcos2, si0, ui3
         vfmsub132pd  xtwo, ur3, ur1
         vfmsub132pd  xtwo, ui3, ui1
         vmovapd      ur0, sr0
         vmovapd      ui0, si0
         vfmadd231pd  tsin1, ui1, sr0
         vfnmadd231pd tsin1, ur1, si0
         vfnmadd132pd tcos1, sr0, ur1
         vfmsub132pd  tcos1, si0, ui1
         vfmsub132pd  xtwo, ur1, ur0
         vfmadd132pd  xtwo, ui1, ui0
         vmovapd      ur2, sr0
         vmovapd      ui2, si0
         vfmadd231pd  tcos1, ui3, sr0
         vfnmadd231pd tcos1, ur3, si0
         vfmadd132pd  tsin1, sr0, ur3
         vfmadd132pd  tsin1, si0, ui3
         vfmsub132pd  xtwo, ur3, ur2
         vfnmadd132pd xtwo, ui3, ui2
         vmovapd      ur0, (%rax)
         vmovapd      ur3, (%rbx)
         vmovapd      ui1, (%rcx)
         vmovapd      ui2, (%rdx)
         pshufb       pbytes, ur2
         pshufb       pbytes, ur1
         pshufb       pbytes, ui3
         pshufb       pbytes, ui0
         vmovupd      ur2, (%rax, %rdi, 8)
         vmovupd      ur1, (%rbx, %rdi, 8)
         vmovupd      ui3, (%rcx, %rdi, 8)
         vmovupd      ui0, (%rdx, %rdi, 8)
         cmp          $8, N2
         jle          done0
         mov          $4, k2
final:   vmovapd      (C, k2, 8), cos1
         vmovapd      (S, k2, 8), sin1
         vmulpd       sin1, sin1, cos2
         vmulpd       cos1, sin1, sin2
         vfmsub231pd  cos1, cos1, cos2
         vfmadd231pd  sin1, cos1, sin2
         lea          (X, k2, 8), %rax
         lea          (%rax, N2, 8), %rbx
         lea          (%rbx, N2, 8), %rcx
         lea          (%rcx, N2, 8), %rdx
         mov          N2, %rdi
         sub          k2, %rdi
         sub          k2, %rdi
         sub          $3, %rdi
         shl          $3, %rdi
         vmovapd      (%rax), er0
         vmovapd      (%rbx), er2
         vmovapd      (%rcx), er1
         vmovapd      (%rdx), er3
         vpermpd      prmt_cntrl, (%rax, %rdi), ei0
         vpermpd      prmt_cntrl, (%rbx, %rdi), ei2
         vpermpd      prmt_cntrl, (%rcx, %rdi), ei1
         vpermpd      prmt_cntrl, (%rdx, %rdi), ei3
         vmovapd      er0, tr0
         vmovapd      ei0, ti0
         vfmadd231pd  sin2, ei2, tr0
         vfnmadd231pd sin2, er2, ti0
         vfnmadd132pd cos2, tr0, er2
         vfnmadd132pd cos2, ti0, ei2
         vfmsub132pd  ytwo, er2, er0
         vfmsub132pd  ytwo, ei2, ei0
         vmovapd      er1, tr0
         vmovapd      ei1, ti0
         vfmadd231pd  sin2, ei3, tr0
         vfnmadd231pd sin2, er3, ti0
         vfnmadd132pd cos2, tr0, er3
         vfnmadd132pd cos2, ti0, ei3
         vfmsub132pd  ytwo, er3, er1
         vfmsub132pd  ytwo, ei3, ei1
         vmovapd      er0, tr0
         vmovapd      ei0, ti0
         vfmadd231pd  sin1, ei1, tr0
         vfnmadd231pd sin1, er1, ti0
         vfnmadd132pd cos1, tr0, er1
         vfmsub132pd  cos1, ti0, ei1
         vfmsub132pd  ytwo, er1, er0
         vfmadd132pd  ytwo, ei1, ei0
         vmovapd      er2, tr0
         vmovapd      ei2, ti0
         vfmadd231pd  cos1, ei3, tr0
         vfnmadd231pd cos1, er3, ti0
         vfmadd132pd  sin1, tr0, er3
         vfmadd132pd  sin1, ti0, ei3
         vfmsub132pd  ytwo, er3, er2
         vfnmadd132pd ytwo, ei3, ei2
         vpermpd      prmt_cntrl, er1, er1
         vpermpd      prmt_cntrl, er2, er2
         vpermpd      prmt_cntrl, ei3, ei3
         vpermpd      prmt_cntrl, ei0, ei0
         vmovapd      er0, (%rax)
         vmovapd      er3, (%rbx)
         vmovapd      ei1, (%rcx)
         vmovapd      ei2, (%rdx)
         vmovupd      er2, (%rax, %rdi)
         vmovupd      er1, (%rbx, %rdi)
         vmovupd      ei3, (%rcx, %rdi)
         vmovupd      ei0, (%rdx, %rdi)
         add          $4, k2
         mov          N2, %rax
         shr          %rax
         cmp          k2, %rax
         jg           final
done0:   pop          %r15
         pop          %r14
         pop          %r13
         pop          %r12
         pop          %rbx
         ret
recurse: bsf          NLO, %rax
         bt           $0, %rax
         jc           radix2
         cmp          $0, k2
         je           keq0
         mov          N2, %rax
         shr          %rax
         cmp          k2, %rax
         je           keqN2o2
         jmp          k2rest
keq0:    xor          ilo, ilo
keq0lp:  lea          (X, ilo, 8), %rax
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
         vmovapd      er2, (%rbx)
         vmovapd      tr2, (%rcx)
         vmovapd      tr3, (%rdx)
         add          $4, ilo
         cmp          ilo, NLO
         jg           keq0lp
         jmp          rcalls
keqN2o2: xor          ilo, ilo
N2o2lp:  lea          (X, ilo, 8), %rax
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
         vmovapd      er2, (%rbx)
         vmovapd      er1, (%rcx)
         vmovapd      er3, (%rdx)
         add          $4, ilo
         cmp          ilo, NLO
         jg           N2o2lp
         jmp          rcalls
k2rest:  mov          N2, %rdx
         dec          %rdx
         mov          X, %rdi
         sub          X0, %rdi
         shr          $5, %rdi
         bsf          NLO, %rcx
         shr          %rcx, %rdi
         mov          %rdi, %rax
         not          %rdi
         and          %rdx, %rdi
         bsf          N2, %rcx
         dec          %rcx
         mov          $1, %rbx
         shl          %rcx, %rbx
L50:     test         %rbx, %rdi
         jz           L100
         xor          %rbx, %rdi
         shr          %rbx
         jmp          L50
L100:    or           %rbx, %rdi
         sub          %rax, %rdi
         jl           rcalls
         imul         NLO, %rdi
         shl          $2, %rdi
         mov          N2, %rax
         shr          %rax
         push         X
         push         k2
         cmp          k2, %rax
         jg           skipneg
         lea          (X, %rdi, 8), X
         neg          %rdi
         mov          N2, %rax
         sub          k2, %rax
         mov          %rax, k2
skipneg: vmovapd      two, ytwo
		 mov          k2, %rax
         imul         NLO, %rax
         imul         $2, %rax, %rbx
         imul         $3, %rax, %rcx
         vbroadcastsd (C, %rax, 8), cos1
         vbroadcastsd (S, %rax, 8), sin1
         vmulpd       sin1, sin1, cos2
         vmulpd       cos1, sin1, sin2
         vfmsub231pd  cos1, cos1, cos2
         vfmadd231pd  sin1, cos1, sin2
         xor          ilo, ilo
k2loop:  lea          (X, ilo, 8), %rax
         lea          (%rax, NLO, 8), %rbx
         lea          (%rbx, NLO, 8), %rcx
         lea          (%rcx, NLO, 8), %rdx
         vmovapd      (%rax), er0
         vmovapd      (%rbx), er1
         vmovapd      (%rcx), er2
         vmovapd      (%rdx), er3
         vmovapd      (%rax, %rdi, 8), ei0
         vmovapd      (%rbx, %rdi, 8), ei1
         vmovapd      (%rcx, %rdi, 8), ei2
         vmovapd      (%rdx, %rdi, 8), ei3
         vmovapd      er0, tr0
         vmovapd      ei0, ti0
         vfmadd231pd  sin2, ei2, tr0
         vfnmadd231pd sin2, er2, ti0
         vfnmadd132pd cos2, tr0, er2
         vfnmadd132pd cos2, ti0, ei2
         vfmsub132pd  ytwo, er2, er0
         vfmsub132pd  ytwo, ei2, ei0
         vmovapd      er1, tr0
         vmovapd      ei1, ti0
         vfmadd231pd  sin2, ei3, tr0
         vfnmadd231pd sin2, er3, ti0
         vfnmadd132pd cos2, tr0, er3
         vfnmadd132pd cos2, ti0, ei3
         vfmsub132pd  ytwo, er3, er1
         vfmsub132pd  ytwo, ei3, ei1
         vmovapd      er0, tr0
         vmovapd      ei0, ti0
         vfmadd231pd  sin1, ei1, tr0
         vfnmadd231pd sin1, er1, ti0
         vfnmadd132pd cos1, tr0, er1
         vfmsub132pd  cos1, ti0, ei1
         vfmsub132pd  ytwo, er1, er0
         vfmadd132pd  ytwo, ei1, ei0
         vmovapd      er2, tr0
         vmovapd      ei2, ti0
         vfmadd231pd  cos1, ei3, tr0
         vfnmadd231pd cos1, er3, ti0
         vfmadd132pd  sin1, tr0, er3
         vfmadd132pd  sin1, ti0, ei3
         vfmsub132pd  ytwo, er3, er2
         vfnmadd132pd ytwo, ei3, ei2
         vmovapd      er0, (%rax)
         vmovapd      ei1, (%rbx)
         vmovapd      er3, (%rcx)
         vmovapd      ei2, (%rdx)
         vmovapd      er2, (%rax, %rdi, 8)
         vmovapd      ei3, (%rbx, %rdi, 8)
         vmovapd      er1, (%rcx, %rdi, 8)
         vmovapd      ei0, (%rdx, %rdi, 8)
         add          $4, ilo
         cmp          ilo, NLO
         jg           k2loop
         pop          k2
         pop          X
rcalls:  cmp          $4, NLO
         jle          done1
         push         X
         push         k2
         bsf          N2, %rcx
         imul         $3, NLO, %rbx
         mov          NLO, %rdx
         shl          %rdx
         neg          %rdx
         mov          $3, %rax
         shl          %rcx, %rax
         or           k2, %rax
         lea          (X, %rbx, 8), X
         push         X
         push         %rax
         mov          $2, %rax
         shl          %rcx, %rax
         or           k2, %rax
         lea          (X, %rdx, 8), X
         push         X
         push         %rax
         mov          $1, %rax
         shl          %rcx, %rax
         or           k2, %rax
         lea          (X, NLO, 8), X
         push         X
         push         %rax
         lea          (X, %rdx, 8), X
         shl          $2, N2
         shr          $2, NLO
         call         recurse
         pop          k2
         pop          X
         call         recurse
         pop          k2
         pop          X
         call         recurse
         pop          k2
         pop          X
         call         recurse
         pop          k2
         pop          X
         shl          $2, NLO
         shr          $2, N2
         jmp          done1
radix2:  xor          ilo, ilo
         shl          NLO
loop_r2: lea          (X, ilo, 8), %rax
         lea          32(X, ilo, 8), %rbx
         lea          64(X, ilo, 8), %rcx
         lea          96(X, ilo, 8), %rdx
         vmovapd      (%rax), %ymm0
         vmovapd      (%rbx), %ymm2
         vmovapd      (%rcx), %ymm4
         vmovapd      (%rdx), %ymm6
         vmovapd      (%rax, NLO, 8), %ymm1
         vmovapd      (%rbx, NLO, 8), %ymm3
         vmovapd      (%rcx, NLO, 8), %ymm5
         vmovapd      (%rdx, NLO, 8), %ymm7
         vaddpd       %ymm1, %ymm0, %ymm8
         vaddpd       %ymm3, %ymm2, %ymm10
         vaddpd       %ymm5, %ymm4, %ymm12
         vaddpd       %ymm7, %ymm6, %ymm14
         vsubpd       %ymm1, %ymm0, %ymm9
         vsubpd       %ymm3, %ymm2, %ymm11
         vsubpd       %ymm5, %ymm4, %ymm13
         vsubpd       %ymm7, %ymm6, %ymm15
         vmovapd      %ymm8, (%rax)
         vmovapd      %ymm10, (%rbx)
         vmovapd      %ymm12, (%rcx)
         vmovapd      %ymm14, (%rdx)
         vmovapd      %ymm9, (%rax, NLO, 8)
         vmovapd      %ymm11, (%rbx, NLO, 8)
         vmovapd      %ymm13, (%rcx, NLO, 8)
         vmovapd      %ymm15, (%rdx, NLO, 8)
         add          $16, ilo
         cmp          ilo, NLO
         jne          loop_r2
         push         X
         push         k2
         bsf          N2, %rcx
         mov          NLO, %rdx
         neg          %rdx
         mov          $1, %rax
         shl          %rcx, %rax
         or           k2, %rax
         lea          (X, NLO, 8), X
         push         X
         push         %rax
         lea          (X, %rdx, 8), X
         shl          N2
         shr          $2, NLO
         call         recurse
         pop          k2
         pop          X
         call         recurse
         pop          k2
         pop          X
         shl          $2, NLO
         shr          N2
done1:   ret





























