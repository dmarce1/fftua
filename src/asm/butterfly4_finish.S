#define  M_SQRT1_2    0.70710678118654752440
#define  N1           $4
#define  SIMD_SIZE    $4
#define  X            %r8
#define  W            %r9
#define  N            %r10
#define  N2           %r11
#define  N2o2         %r12
#define  _8N2o2       %r13
#define  k2max        %r14
#define  k2           %r15
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
#define  sr0          %xmm4
#define  sr1          %xmm5
#define  sr2          %xmm6
#define  sr3          %xmm7
#define  ui0          %xmm8
#define  ui1          %xmm9
#define  ui2          %xmm10
#define  ui3          %xmm11
#define  si0          %xmm12
#define  si1          %xmm13
#define  si2          %xmm14
#define  si3          %xmm15
#define  prmt_cntrl   $27

         .global      butterfly4_finish

         .data

         .align       32
dI1:     .double      0
         .double      2
         .double      4
         .double      6
dI2:     .double      0
         .double      4
         .double      8
         .double      12
dI3:     .double      0
         .double      6
         .double      12
         .double      18
ones:    .quad        0xffffffffffffffff
         .quad        0xffffffffffffffff
         .quad        0xffffffffffffffff
         .quad        0xffffffffffffffff
tw45:    .double      M_SQRT1_2
none:    .double      -1.0

         .text

butterfly4_finish:
         push         %rbx
         push         %r15
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
         movq         (%rax), ur0
         movq         (%rbx), ur1
         movq         (%rcx), ur2
         movq         (%rdx), ur3
         vaddsd       ur2, ur0, sr0
         vsubsd       ur2, ur0, sr2
         vaddsd       ur3, ur1, sr1
         vsubsd       ur3, ur1, sr3
         vaddsd       sr1, sr0, ur0
         vsubsd       sr1, sr0, ur2
         movq         ur0, (%rax)
         movq         sr2, (%rbx)
         movq         sr3, (%rcx)
         movq         ur2, (%rdx)
         cmp          $1, N2
         jle          done
         movq         (%rax, _8N2o2), ur0
         movq         (%rbx, _8N2o2), ur1
         movq         (%rcx, _8N2o2), ur2
         movq         (%rdx, _8N2o2), ur3
         vaddsd       ur3, ur1, sr0
         vsubsd       ur3, ur1, sr2
         vmulsd       tw45, sr2, sr1
         vmulsd       tw45, sr0, sr3
         movq         ur0, sr0
         movq         ur2, sr2
         vaddsd       sr1, sr0, ur0
         vaddsd       sr3, sr2, ur3
         vsubsd       sr1, sr0, ur1
         vsubsd       sr3, sr2, ur2
         vmulsd       none, ur3, ur3
         movq         ur0, (%rax, _8N2o2)
         movq         ur1, (%rbx, _8N2o2)
         movq         ur2, (%rcx, _8N2o2)
         movq         ur3, (%rdx, _8N2o2)
         cmp          $2, N2
         jle          done
         movq         $8, k2max
         cmp          $8, N2
         cmovl        N2, k2max
         shr          k2max
         mov          $1, k2
k2scalar:
         mov          k2, %rax
         shl          %rax
         mov          %rax, %rbx
         mov          %rax, %rcx
         add          %rax, %rbx
         add          %rbx, %rcx
         movq         0(W, %rax, 8), sr1
         movq         0(W, %rbx, 8), sr2
         movq         0(W, %rcx, 8), sr3
         movq         8(W, %rax, 8), si1
         movq         8(W, %rbx, 8), si2
         movq         8(W, %rcx, 8), si3
         movq         N, %rdi
         sub          k2, %rdi
         sub          k2, %rdi
         lea          (X, k2, 8), %rax
         lea          (%rax, N2, 8), %rbx
         lea          (%rbx, N2, 8), %rcx
         lea          (%rcx, N2, 8), %rdx
         movq         (%rax), ur0
         movq         (%rbx), ur1
         movq         (%rcx), ur2
         movq         (%rdx), ur3
         movq         (%rax, %rdi, 8), ui0
         movq         (%rbx, %rdi, 8), ui1
         movq         (%rcx, %rdi, 8), ui2
         movq         (%rdx, %rdi, 8), ui3
         vmulsd       si1, ui1, sr0
         vmulsd       sr1, ui1, si0
         vfmsub231sd  sr1, ur1, sr0
         vfmadd231sd  si1, ur1, si0
         movq         si0, ui1
         movq         sr0, ur1
         vmulsd       si2, ui2, sr0
         vmulsd       sr2, ui2, si0
         vfmsub231sd  sr2, ur2, sr0
         vfmadd231sd  si2, ur2, si0
         movq         si0, ui2
         movq         sr0, ur2
         vmulsd       si3, ui3, sr0
         vmulsd       sr3, ui3, si0
         vfmsub231sd  sr3, ur3, sr0
         vfmadd231sd  si3, ur3, si0
         movq         si0, ui3
         movq         sr0, ur3
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
         movq         (%rax), ur0
         movq         (%rbx), ur1
         movq         (%rcx), ui2
         movq         (%rdx), ui3
         movq         (%rax, %rdi, 8), ur3
         movq         (%rbx, %rdi, 8), ur2
         movq         (%rcx, %rdi, 8), ui1
         movq         (%rdx, %rdi, 8), ui0
         inc          k2
         cmp          k2, k2max
         jne          k2scalar
         cmp          $4, N2
         jle          done
         mov          $4, k2
k2simd:
         mov          k2, %rax
         shl          %rax
         mov          %rax, %rbx
         mov          %rax, %rcx
         add          %rax, %rbx
         add          %rbx, %rcx
         lea          (W, %rax, 8), %rax
         lea          (W, %rbx, 8), %rbx
         lea          (W, %rcx, 8), %rcx
         movd         dI1, %xmm1
         movd         dI2, %xmm2
         movd         dI3, %xmm3
         vmovdqa      ones, %ymm4
         vmovdqa      %ymm4, %ymm5
         vmovdqa      %ymm4, %ymm6
         vmovdqa      %ymm4, %ymm7
         vmovdqa      %ymm4, %ymm8
         vmovdqa      %ymm4, %ymm9
         vgatherdpd   %ymm4, 0(%rax, %xmm1, 8), tr1
         vgatherdpd   %ymm5, 8(%rax, %xmm1, 8), ti1
         vgatherdpd   %ymm6, 0(%rbx, %xmm2, 8), tr2
         vgatherdpd   %ymm7, 8(%rbx, %xmm2, 8), ti2
         vgatherdpd   %ymm8, 0(%rcx, %xmm3, 8), tr3
         vgatherdpd   %ymm9, 8(%rcx, %xmm3, 8), ti3
         mov          N2, %rdi
         sub          k2, %rdi
         sub          k2, %rdi
         sub          $24, %rdi
         lea          (W, k2, 8), %rax
         lea          (%rax, N2, 8), %rbx
         lea          (%rbx, N2, 8), %rcx
         lea          (%rcx, N2, 8), %rdx
         vmovapd      (%rax), er0
         vmovapd      (%rbx), er1
         vmovapd      (%rcx), er2
         vmovapd      (%rdx), er3
         vmovupd      (%rax, %rdi, 8), ei0
         vmovupd      (%rbx, %rdi, 8), ei1
         vmovupd      (%rcx, %rdi, 8), ei2
         vmovupd      (%rdx, %rdi, 8), ei3
         vpermpd      prmt_cntrl, ei0, ei0
         vpermpd      prmt_cntrl, ei1, ei1
         vpermpd      prmt_cntrl, ei2, ei2
         vpermpd      prmt_cntrl, ei3, ei3
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
         vmovupd      er3, (%rax, %rdi, 8)
         vmovupd      er2, (%rbx, %rdi, 8)
         vmovupd      ei1, (%rcx, %rdi, 8)
         vmovupd      ei0, (%rdx, %rdi, 8)
         add          SIMD_SIZE, k2
         cmp          k2, N2o2
         jne          k2simd
done:
         pop          %r12
         pop          %r13
         pop          %r14
         pop          %r15
         pop          %rbx
         ret
