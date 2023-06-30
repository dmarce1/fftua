#define       X              %r8
#define       C              %r9
#define       N              %r10
#define       N2             %r11
#define       k2             %r12
#define       er0            %ymm0
#define       er1            %ymm1
#define       er2            %ymm2
#define       er3            %ymm3
#define       ei0            %ymm4
#define       ei1            %ymm5
#define       ei2            %ymm6
#define       ei3            %ymm7
#define       tr0            %ymm8
#define       ti0            %ymm9
#define       tr1            %ymm10
#define       ti1            %ymm11
#define       tr2            %ymm12
#define       ti2            %ymm13
#define       tr3            %ymm14
#define       ti3            %ymm15
#define       ur0            %xmm0
#define       ui0            %xmm1
#define       ur1            %xmm2
#define       ui1            %xmm3
#define       ur2            %xmm4
#define       ui2            %xmm5
#define       ur3            %xmm6
#define       ui3            %xmm7
#define       sr0            %xmm8
#define       si0            %xmm9
#define       sr1            %xmm10
#define       sr2            %xmm11
#define       si1            %xmm12
#define       si2            %xmm13
#define       sr3            %xmm14
#define       si3            %xmm15

              .global        dit_rn_final

              .data

              .align         32
laneswap16:   .byte          8
              .byte          9
              .byte          10
              .byte          11
              .byte          12
              .byte          13
              .byte          14
              .byte          15
              .byte          0
              .byte          1
              .byte          2
              .byte          3
              .byte          4
              .byte          5
              .byte          6
              .byte          7

              .text

dit_rn_final: push           %rbx
              push           %r12
              mov            %rdi, X
              mov            %rsi, C
              mov            %rdx, N
              mov            N, N2
              shr            $2, N2
              lea            (X), %rax
              lea            (%rax, N2, 8), %rbx
              lea            (%rbx, N2, 8), %rcx
              lea            (%rcx, N2, 8), %rdx
              movsd          (%rax), ur0
              movsd          (%rbx), ur2
              movsd          (%rcx), ur1
              movsd          (%rdx), ur3
              vaddsd         ur2, ur0, sr0
              vsubsd         ur2, ur0, sr2
              vaddsd         ur3, ur1, sr1
              vsubsd         ur1, ur3, sr3
              vaddsd         sr1, sr0, ur0
              vsubsd         sr1, sr0, ur2
              movsd          ur0, (%rax)
              movsd          sr2, (%rbx)
              movsd          ur2, (%rcx)
              movsd          sr3, (%rdx)
              cmp            $1, N2
              jle            done0
              movsd          (%rax, N), ur0
              movsd          (%rbx, N), ur2
              movsd          (%rcx, N), ur1
              movsd          (%rdx, N), ur3
              vaddsd         ur3, ur1, sr0
              vsubsd         ur3, ur1, sr2
              vmulsd         tw45, sr2, sr1
              vmulsd         tw45, sr0, sr3
              movsd          ur0, sr0
              movsd          ur2, sr2
              vaddsd         sr1, sr0, ur0
              vaddsd         sr3, sr2, ur3
              vsubsd         sr1, sr0, ur1
              vsubsd         sr3, sr2, ur2
              vmulsd         none, ur3, ur3
              movsd          ur0, (%rax, N)
              movsd          ur1, (%rbx, N)
              movsd          ur2, (%rcx, N)
              movsd          ur3, (%rdx, N)
              cmp            $2, N2
              jle            done0
              movsd          8(C), sr1
              movsd          8(C, N, 2), si1
              vmulsd         si1, si1, sr2
              vmulsd         sr1, si1, si2
              vfmsub231sd    sr1, sr1, sr2
              vfmadd231sd    si1, sr1, si2
              mov            N2, %rdi
              sub            $2, %rdi
              lea            8(X), %rax
              lea            (%rax, N2, 8), %rbx
              lea            (%rbx, N2, 8), %rcx
              lea            (%rcx, N2, 8), %rdx
              movsd          (%rax), ur0
              movsd          (%rbx), ur2
              movsd          (%rcx), ur1
              movsd          (%rdx), ur3
              movsd          (%rax, %rdi, 8), ui0
              movsd          (%rbx, %rdi, 8), ui2
              movsd          (%rcx, %rdi, 8), ui1
              movsd          (%rdx, %rdi, 8), ui3
              movsd          ur0, sr0
              movsd          ui0, si0
              vfmadd231sd    si2, ui2, sr0
              vfnmadd231sd   si2, ur2, si0
              vfnmadd132sd   sr2, sr0, ur2
              vfnmadd132sd   sr2, si0, ui2
              vfmsub132sd    two, ur2, ur0
              vfmsub132sd    two, ui2, ui0
              movsd          ur1, sr0
              movsd          ui1, si0
              vfmadd231sd    si2, ui3, sr0
              vfnmadd231sd   si2, ur3, si0
              vfnmadd132sd   sr2, sr0, ur3
              vfnmadd132sd   sr2, si0, ui3
              vfmsub132sd    two, ur3, ur1
              vfmsub132sd    two, ui3, ui1
              movsd          ur0, sr0
              movsd          ui0, si0
              vfmadd231sd    si1, ui1, sr0
              vfnmadd231sd   si1, ur1, si0
              vfnmadd132sd   sr1, sr0, ur1
              vfmsub132sd    sr1, si0, ui1
              vfmsub132sd    two, ur1, ur0
              vfmadd132sd    two, ui1, ui0
              movsd          ur2, sr0
              movsd          ui2, si0
              vfmadd231sd    sr1, ui3, sr0
              vfnmadd231sd   sr1, ur3, si0
              vfmadd132sd    si1, sr0, ur3
              vfmadd132sd    si1, si0, ui3
              vfmsub132sd    two, ur3, ur2
              vfnmadd132sd   two, ui3, ui2
              movsd          ur0, (%rax)
              movsd          ur3, (%rbx)
              movsd          ui1, (%rcx)
              movsd          ui2, (%rdx)
              movsd          ur2, (%rax, %rdi, 8)
              movsd          ur1, (%rbx, %rdi, 8)
              movsd          ui3, (%rcx, %rdi, 8)
              movsd          ui0, (%rdx, %rdi, 8)
              cmp            $4, N2
              jle            done0
              movlpd         16(C), sr1
              movlpd         16(C, N, 2), si1
              vmovhpd        24(C), sr1, sr1
              vmovhpd        24(C, N, 2), si1, si1
              vmulpd         si1, si1, sr2
              vmulpd         sr1, si1, si2
              vfmsub231pd    sr1, sr1, sr2
              vfmadd231pd    si1, sr1, si2
              mov            N2, %rdi
              sub            $5, %rdi
              lea            16(X), %rax
              lea            (%rax, N2, 8), %rbx
              lea            (%rbx, N2, 8), %rcx
              lea            (%rcx, N2, 8), %rdx
              vmovapd        (%rax), ur0
              vmovapd        (%rbx), ur2
              vmovapd        (%rcx), ur1
              vmovapd        (%rdx), ur3
              vmovupd        (%rax, %rdi, 8), ui0
              vmovupd        (%rbx, %rdi, 8), ui2
              vmovupd        (%rcx, %rdi, 8), ui1
              vmovupd        (%rdx, %rdi, 8), ui3
              pshufb         laneswap16, ui0
              pshufb         laneswap16, ui1
              pshufb         laneswap16, ui2
              pshufb         laneswap16, ui3
              vmovapd        ur0, sr0
              vmovapd        ui0, si0
              vfmadd231pd    si2, ui2, sr0
              vfnmadd231pd   si2, ur2, si0
              vfnmadd132pd   sr2, sr0, ur2
              vfnmadd132pd   sr2, si0, ui2
              vfmsub132pd    two, ur2, ur0
              vfmsub132pd    two, ui2, ui0
              vmovapd        ur1, sr0
              vmovapd        ui1, si0
              vfmadd231pd    si2, ui3, sr0
              vfnmadd231pd   si2, ur3, si0
              vfnmadd132pd   sr2, sr0, ur3
              vfnmadd132pd   sr2, si0, ui3
              vfmsub132pd    two, ur3, ur1
              vfmsub132pd    two, ui3, ui1
              vmovapd        ur0, sr0
              vmovapd        ui0, si0
              vfmadd231pd    si1, ui1, sr0
              vfnmadd231pd   si1, ur1, si0
              vfnmadd132pd   sr1, sr0, ur1
              vfmsub132pd    sr1, si0, ui1
              vfmsub132pd    two, ur1, ur0
              vfmadd132pd    two, ui1, ui0
              vmovapd        ur2, sr0
              vmovapd        ui2, si0
              vfmadd231pd    sr1, ui3, sr0
              vfnmadd231pd   sr1, ur3, si0
              vfmadd132pd    si1, sr0, ur3
              vfmadd132pd    si1, si0, ui3
              vfmsub132pd    two, ur3, ur2
              vfnmadd132pd   two, ui3, ui2
              vmovapd        ur0, (%rax)
              vmovapd        ur3, (%rbx)
              vmovapd        ui1, (%rcx)
              vmovapd        ui2, (%rdx)
              pshufb         laneswap16, ur2
              pshufb         laneswap16, ur1
              pshufb         laneswap16, ui3
              pshufb         laneswap16, ui0
              vmovupd        ur2, (%rax, %rdi, 8)
              vmovupd        ur1, (%rbx, %rdi, 8)
              vmovupd        ui3, (%rcx, %rdi, 8)
              vmovupd        ui0, (%rdx, %rdi, 8)
              cmp            $8, N2
              jle            done0
              mov            $4, k2
final:        lea            (C, k2, 8), %rax
              vmovapd        (%rax), tr1
              vmovapd        (%rax, N, 2), ti1
              vmulpd         ti1, ti1, tr2
              vmulpd         tr1, ti1, ti2
              vfmsub231pd    tr1, tr1, tr2
              vfmadd231pd    ti1, tr1, ti2
              lea            (X, k2, 8), %rax
              lea            (%rax, N2, 8), %rbx
              lea            (%rbx, N2, 8), %rcx
              lea            (%rcx, N2, 8), %rdx
              mov            N2, %rdi
              vmovapd        (%rax), er0
              sub            k2, %rdi
              vmovapd        (%rbx), er2
              sub            k2, %rdi
              vmovapd        (%rcx), er1
              sub            $3, %rdi
              vmovapd        (%rdx), er3
              shl            $3, %rdi
              vpermpd        $27, (%rax, %rdi), ei0
              vpermpd        $27, (%rbx, %rdi), ei2
              vpermpd        $27, (%rcx, %rdi), ei1
              vpermpd        $27, (%rdx, %rdi), ei3
              vmovapd        er0, tr0
              vmovapd        ei0, ti0
              vfmadd231pd    ti2, ei2, tr0
              vfnmadd231pd   ti2, er2, ti0
              vfnmadd132pd   tr2, tr0, er2
              vfnmadd132pd   tr2, ti0, ei2
              vfmsub132pd    two, er2, er0
              vfmsub132pd    two, ei2, ei0
              vmovapd        er1, tr0
              vmovapd        ei1, ti0
              vfmadd231pd    ti2, ei3, tr0
              vfnmadd231pd   ti2, er3, ti0
              vfnmadd132pd   tr2, tr0, er3
              vfnmadd132pd   tr2, ti0, ei3
              vfmsub132pd    two, er3, er1
              vfmsub132pd    two, ei3, ei1
              vmovapd        er0, tr0
              vmovapd        ei0, ti0
              vfmadd231pd    ti1, ei1, tr0
              vfnmadd231pd   ti1, er1, ti0
              vfnmadd132pd   tr1, tr0, er1
              vfmsub132pd    tr1, ti0, ei1
              vfmsub132pd    two, er1, er0
              vfmadd132pd    two, ei1, ei0
              vmovapd        er2, tr0
              vmovapd        ei2, ti0
              vfmadd231pd    tr1, ei3, tr0
              vfnmadd231pd   tr1, er3, ti0
              vfmadd132pd    ti1, tr0, er3
              vfmadd132pd    ti1, ti0, ei3
              vfmsub132pd    two, er3, er2
              vfnmadd132pd   two, ei3, ei2
              vpermpd        $27, er1, er1
              vpermpd        $27, er2, er2
              vpermpd        $27, ei3, ei3
              vpermpd        $27, ei0, ei0
              vmovapd        er0, (%rax)
              vmovapd        er3, (%rbx)
              vmovapd        ei1, (%rcx)
              vmovapd        ei2, (%rdx)
              vmovupd        er2, (%rax, %rdi)
              vmovupd        er1, (%rbx, %rdi)
              vmovupd        ei3, (%rcx, %rdi)
              vmovupd        ei0, (%rdx, %rdi)
              add            $4, k2
              mov            N2, %rax
              shr            %rax
              cmp            k2, %rax
              jg             final
done0:        pop            %r12
              pop            %rbx
              ret
