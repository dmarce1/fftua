#define  M_SQRT1_2    0.70710678118654752440
#define  N1           $4
#define  SIMD_SIZE    $4
#define  X            %r8
#define  W            %r9
#define  N            %r10
#define  N2           %r11
#define  N2o2         %r12
#define  _8N2o2       %r13
#define  k2           %r14
#define  er0          %ymm0
#define  er1          %ymm1
#define  er2          %ymm2
#define  er3          %ymm3
#define  tr0          %ymm4
#define  ei0          %ymm5
#define  ei1          %ymm6
#define  ei2          %ymm7
#define  ei3          %ymm8
#define  ti0          %ymm9
#define  tr1          %ymm10
#define  tr2          %ymm11
#define  tr3          %ymm12
#define  ti1          %ymm13
#define  ti2          %ymm14
#define  ti3          %ymm15
#define  ur0          %xmm0
#define  ur1          %xmm1
#define  ur2          %xmm2
#define  ur3          %xmm3
#define  ui0          %xmm4
#define  ui1          %xmm5
#define  ui2          %xmm6
#define  ui3          %xmm7
#define  sr0          %xmm8
#define  sr1          %xmm9
#define  sr2          %xmm10
#define  sr3          %xmm11
#define  si0          %xmm12
#define  si1          %xmm13
#define  si2          %xmm14
#define  si3          %xmm15
#define  prmt_cntrl   $27

         .global      butterfly4_finish

         .data
         .align       32
ones:    .quad        0xffffffffffffffff
         .quad        0xffffffffffffffff
         .quad        0xffffffffffffffff
         .quad        0xffffffffffffffff
dI1:     .quad        0
         .quad        2
         .quad        4
         .quad        6
dI2:     .quad        0
         .quad        4
         .quad        8
         .quad        12
dI3:     .quad        0
         .quad        6
         .quad        12
         .quad        18
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
tw45:    .double      M_SQRT1_2
none:    .double      -1.0

         .text

butterfly4_finish:
         push         %rbx
         push         %r14
         push         %r13
         push         %r12
         mov          %rdi, X
         mov          %rsi, W
         mov          %rdx, N
         shr          $2, %rdx
         mov          %rdx, N2
         shl          $2, %rdx
         mov          %rdx, _8N2o2
         shr          $3, %rdx
         mov          %rdx, N2o2
         lea          (X), %rax
         lea          (%rax, N2, 8), %rbx
         lea          (%rbx, N2, 8), %rcx
         lea          (%rcx, N2, 8), %rdx
         movsd        (%rax), ur0
         movsd        (%rbx), ur1
         movsd        (%rcx), ur2
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
         jle          done
         movsd        (%rax, _8N2o2), ur0
         movsd        (%rbx, _8N2o2), ur1
         movsd        (%rcx, _8N2o2), ur2
         movsd        (%rdx, _8N2o2), ur3
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
         movsd        ur0, (%rax, _8N2o2)
         movsd        ur1, (%rbx, _8N2o2)
         movsd        ur2, (%rcx, _8N2o2)
         movsd        ur3, (%rdx, _8N2o2)
         cmp          $2, N2
         jle          done
         movsd        16(W), sr1
         movsd        32(W), sr2
         movsd        48(W), sr3
         movsd        24(W), si1
         movsd        40(W), si2
         movsd        56(W), si3
         mov          N2, %rdi
         sub          $2, %rdi
         lea          8(X), %rax
         lea          (%rax, N2, 8), %rbx
         lea          (%rbx, N2, 8), %rcx
         lea          (%rcx, N2, 8), %rdx
         movsd        (%rax), ur0
         movsd        (%rbx), ur1
         movsd        (%rcx), ur2
         movsd        (%rdx), ur3
         movsd        (%rax, %rdi, 8), ui0
         movsd        (%rbx, %rdi, 8), ui1
         movsd        (%rcx, %rdi, 8), ui2
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
         jle          done
         mov          $2, k2
         movlpd       32(W), sr1
         movlpd       64(W), sr2
         movlpd       96(W), sr3
         movlpd       40(W), si1
         movlpd       72(W), si2
         movlpd       104(W), si3
         vmovhpd      48(W), sr1, sr1
         vmovhpd      96(W), sr2, sr2
         vmovhpd      144(W), sr3, sr3
         vmovhpd      56(W), si1, si1
         vmovhpd      104(W), si2, si2
         vmovhpd      152(W), si3, si3
         mov          N2, %rdi
         sub          $5, %rdi
         lea          16(X), %rax
         lea          (%rax, N2, 8), %rbx
         lea          (%rbx, N2, 8), %rcx
         lea          (%rcx, N2, 8), %rdx
         vmovapd      (%rax), ur0
         vmovapd      (%rbx), ur1
         vmovapd      (%rcx), ur2
         vmovapd      (%rdx), ur3
         vmovupd      (%rax, %rdi, 8), ui0
         vmovupd      (%rbx, %rdi, 8), ui1
         vmovupd      (%rcx, %rdi, 8), ui2
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
         jle          done
         mov          $4, k2
k2simd:
         mov          k2, %rax
         shl          %rax
         lea          (W, %rax, 8), %rax
         vmovdqa      dI1, %ymm1
         vmovdqa      ones, %ymm4
         vmovdqa      %ymm4, %ymm5
         vgatherqpd   %ymm4, 0(%rax, %ymm1, 8), tr1
         vgatherqpd   %ymm5, 8(%rax, %ymm1, 8), ti1
         vmulpd       ti1, ti1, tr2
         vmulpd       tr1, ti1, ti2
         vfmsub231pd  tr1, tr1, tr2
         vfmadd231pd  ti1, tr1, ti2
         vmulpd       ti1, ti2, tr3
         vmulpd       tr1, ti2, ti3
         vfmsub231pd  tr1, tr2, tr3
         vfmadd231pd  ti1, tr2, ti3
         mov          N2, %rdi
         sub          k2, %rdi
         sub          k2, %rdi
         sub          $3, %rdi
         lea          (X, k2, 8), %rax
         lea          (%rax, N2, 8), %rbx
         lea          (%rbx, N2, 8), %rcx
         lea          (%rcx, N2, 8), %rdx
         vmovapd      (%rax), er0
         vmovapd      (%rbx), er1
         vmovapd      (%rcx), er2
         vmovapd      (%rdx), er3
         vpermpd      prmt_cntrl, (%rax, %rdi, 8), ei0
         vpermpd      prmt_cntrl, (%rbx, %rdi, 8), ei1
         vpermpd      prmt_cntrl, (%rcx, %rdi, 8), ei2
         vpermpd      prmt_cntrl, (%rdx, %rdi, 8), ei3
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
         vmovapd      ei3, (%rdx)
         vmovapd      ei2, (%rcx)
         vmovupd      er2, (%rbx, %rdi, 8)
         vmovupd      er3, (%rax, %rdi, 8)
         vmovupd      ei1, (%rcx, %rdi, 8)
         vmovupd      ei0, (%rdx, %rdi, 8)
         add          SIMD_SIZE, k2
         cmp          k2, N2o2
         jne          k2simd
done:
         pop          %r12
         pop          %r13
         pop          %r14
         pop          %rbx
         ret
