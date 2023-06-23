#define  M_SQRT1_2    0.70710678118654752440
#define  N1           $4
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
#define  tr1          %ymm9
#define  tr2          %ymm10
#define  tr3          %ymm11
#define  ti0          %ymm12
#define  ti1          %ymm13
#define  ti2          %ymm14
#define  ti3          %ymm15
#define  ur0          %xmm0
#define  ui0          %xmm1
#define  ur1          %xmm2
#define  ui1          %xmm3
#define  ur2          %xmm4
#define  ui2          %xmm5
#define  ur3          %xmm6
#define  ui3          %xmm7
#define  sr0          %xmm8
#define  sr1          %xmm9
#define  sr2          %xmm10
#define  sr3          %xmm11
#define  si0          %xmm12
#define  si1          %xmm13
#define  si2          %xmm14
#define  si3          %xmm15
#define  cos1         -8(%rbp)
#define  cos2         -16(%rbp)
#define  cos3         -24(%rbp)
#define  sin1         -32(%rbp)
#define  sin2         -40(%rbp)
#define  sin3         -48(%rbp)
#define  X0           -56(%rbp)
#define  prmt_cntrl   $27
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
trash:   .double      -99.0
         .double      -99.0
         .double      -99.0
         .double      -99.0
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
         enter        STACK_SIZE, $0
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
         movsd        sr2, (%rbx)
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
         movsd        ur2, sr2
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
         movsd        8(C), sr1
         movsd        16(C), sr2
         movsd        24(C), sr3
         movsd        8(S), si1
         movsd        16(S), si2
         movsd        24(S), si3
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
         vmulsd       si1, ui1, sr0
         vmulsd       sr1, ui1, si0
         vfmsub231sd  sr1, ur1, sr0
         vfmadd231sd  si1, ur1, si0
         movsd        si0, ui1
         movsd        sr0, ur1
         vmulsd       si2, ui2, sr0
         vmulsd       sr2, ui2, si0
         vfmsub231sd  sr2, ur2, sr0
         vfmadd231sd  si2, ur2, si0
         movsd        si0, ui2
         movsd        sr0, ur2
         vmulsd       si3, ui3, sr0
         vmulsd       sr3, ui3, si0
         vfmsub231sd  sr3, ur3, sr0
         vfmadd231sd  si3, ur3, si0
         movsd        si0, ui3
         movsd        sr0, ur3
         vaddsd       ur2, ur0, sr0
         vaddsd       ui2, ui0, si0
         vsubsd       ur2, ur0, sr2
         vsubsd       ui2, ui0, si2
         vaddsd       ur3, ur1, sr1
         vaddsd       ui3, ui1, si1
         vsubsd       ur1, ur3, sr3
         vsubsd       ui3, ui1, si3
         vaddsd       sr1, sr0, ur0
         vaddsd       si1, si0, ui0
         vaddsd       si3, sr2, ur1
         vaddsd       sr3, si2, ui1
         vsubsd       sr1, sr0, ur2
         vsubsd       si0, si1, ui2
         vsubsd       si3, sr2, ur3
         vsubsd       si2, sr3, ui3
         movsd        ur0, (%rax)
         movsd        ur1, (%rbx)
         movsd        ui2, (%rcx)
         movsd        ui3, (%rdx)
         movsd        ur3, (%rax, %rdi, 8)
         movsd        ur2, (%rbx, %rdi, 8)
         movsd        ui1, (%rcx, %rdi, 8)
         movsd        ui0, (%rdx, %rdi, 8)
         cmp          $4, N2
         jle          done0
         movlpd       16(C), sr1
         movlpd       32(C), sr2
         movlpd       48(C), sr3
         movlpd       16(S), si1
         movlpd       32(S), si2
         movlpd       48(S), si3
         vmovhpd      24(C), sr1, sr1
         vmovhpd      48(C), sr2, sr2
         vmovhpd      72(C), sr3, sr3
         vmovhpd      24(S), si1, si1
         vmovhpd      48(S), si2, si2
         vmovhpd      72(S), si3, si3
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
         vmulpd       si1, ui1, sr0
         vmulpd       sr1, ui1, si0
         vfmsub231pd  sr1, ur1, sr0
         vfmadd231pd  si1, ur1, si0
         vmovapd      si0, ui1
         vmovapd      sr0, ur1
         vmulpd       si2, ui2, sr0
         vmulpd       sr2, ui2, si0
         vfmsub231pd  sr2, ur2, sr0
         vfmadd231pd  si2, ur2, si0
         vmovapd      si0, ui2
         vmovapd      sr0, ur2
         vmulpd       si3, ui3, sr0
         vmulpd       sr3, ui3, si0
         vfmsub231pd  sr3, ur3, sr0
         vfmadd231pd  si3, ur3, si0
         vmovapd      si0, ui3
         vmovapd      sr0, ur3
         vaddpd       ur2, ur0, sr0
         vaddpd       ui2, ui0, si0
         vsubpd       ur2, ur0, sr2
         vsubpd       ui2, ui0, si2
         vaddpd       ur3, ur1, sr1
         vaddpd       ui3, ui1, si1
         vsubpd       ur1, ur3, sr3
         vsubpd       ui3, ui1, si3
         vaddpd       sr1, sr0, ur0
         vaddpd       si1, si0, ui0
         vaddpd       si3, sr2, ur1
         vaddpd       sr3, si2, ui1
         vsubpd       sr1, sr0, ur2
         vsubpd       si0, si1, ui2
         vsubpd       si3, sr2, ur3
         vsubpd       si2, sr3, ui3
         vmovapd      ur0, (%rax)
         vmovapd      ur1, (%rbx)
         vmovapd      ui2, (%rcx)
         vmovapd      ui3, (%rdx)
         pshufb       pbytes, ur3
         pshufb       pbytes, ur2
         pshufb       pbytes, ui1
         pshufb       pbytes, ui0
         vmovupd      ur3, (%rax, %rdi, 8)
         vmovupd      ur2, (%rbx, %rdi, 8)
         vmovupd      ui1, (%rcx, %rdi, 8)
         vmovupd      ui0, (%rdx, %rdi, 8)
         cmp          $8, N2
         jle          done0
         mov          $4, k2
final:   vmovapd      (C, k2, 8), tr1
         vmovapd      (S, k2, 8), ti1
         vmulpd       ti1, ti1, tr2
         vmulpd       tr1, ti1, ti2
         vfmsub231pd  tr1, tr1, tr2
         vfmadd231pd  ti1, tr1, ti2
         vmulpd       ti1, ti2, tr3
         vmulpd       tr1, ti2, ti3
         vfmsub231pd  tr1, tr2, tr3
         vfmadd231pd  ti1, tr2, ti3
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
         vpermpd      prmt_cntrl, er3, er3
         vpermpd      prmt_cntrl, er2, er2
         vpermpd      prmt_cntrl, ei1, ei1
         vpermpd      prmt_cntrl, ei0, ei0
         vmovapd      er0, (%rax)
         vmovapd      er1, (%rbx)
         vmovapd      ei2, (%rcx)
         vmovapd      ei3, (%rdx)
         vmovupd      er3, (%rax, %rdi)
         vmovupd      er2, (%rbx, %rdi)
         vmovupd      ei1, (%rcx, %rdi)
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
         leave
         ret
recurse: cmp          $1, NLO
         je           done1
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
skipneg: mov          k2, %rax
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
         xor          ilo, ilo
k2loop:  vbroadcastsd cos1, tr1
         vbroadcastsd cos2, tr2
         vbroadcastsd cos3, tr3
         vbroadcastsd sin1, ti1
         vbroadcastsd sin2, ti2
         vbroadcastsd sin3, ti3
         lea          (X, ilo, 8), %rax
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
  //       er0 er1 er3 er2 ei2 ei3 ei1 ei0
  //       er0 ei2 er3 ei1 er1 ei3 er2 ei0
   //      er0 ei0 er1 ei1 er2 ei2 er3 ei3
//         vmovapd      trash, er3
         vmovapd      er0, (%rax)
         vmovapd      ei2, (%rbx)
         vmovapd      er1, (%rcx)
         vmovapd      ei3, (%rdx)
         vmovapd      er3, (%rax, %rdi, 8)
         vmovapd      ei1, (%rbx, %rdi, 8)
         vmovapd      er2, (%rcx, %rdi, 8)
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
done1:   ret





























