#include      "common.h"


              .global        fft_recursive


              .text

fft_recursive:push           %rbx
              push           %r12
              push           %r13
              push           %r14
              push           %r15
              mov            %rdi, X
              mov            %rsi, C
              mov            %rdx, S
              mov            %rcx, N0
              mov            N0, NLO
              mov            N0, N
              shr            $2, NLO
              mov            $1, N2
              xor            k2, k2
              mov            X, X0

              call           fft_dispatch

              mov            X, %rdi
              mov            N0, %rsi
              push           %r8
              push           %r9
              push           %r10
              push           %r11
              call           fft_scramble
              pop            %r11
              pop            %r10
              pop            %r9
              pop            %r8
              mov            N0, N2
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
              movsd          (%rax, N0), ur0
              movsd          (%rbx, N0), ur2
              movsd          (%rcx, N0), ur1
              movsd          (%rdx, N0), ur3
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
              movsd          ur0, (%rax, N0)
              movsd          ur1, (%rbx, N0)
              movsd          ur2, (%rcx, N0)
              movsd          ur3, (%rdx, N0)
              cmp            $2, N2
              jle            done0
              movsd          8(C), sr1
              movsd          8(S), si1
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
              movlpd         16(S), si1
              vmovhpd        24(C), sr1, sr1
              vmovhpd        24(S), si1, si1
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
final:        vmovapd        (C, k2, 8), tr1
              vmovapd        (S, k2, 8), ti1
              vmulpd         ti1, ti1, tr2
              vmulpd         tr1, ti1, ti2
              vfmsub231pd    tr1, tr1, tr2
              vfmadd231pd    ti1, tr1, ti2
              lea            (X, k2, 8), %rax
              lea            (%rax, N2, 8), %rbx
              lea            (%rbx, N2, 8), %rcx
              lea            (%rcx, N2, 8), %rdx
              mov            N2, %rdi
              sub            k2, %rdi
              sub            k2, %rdi
              sub            $3, %rdi
              shl            $3, %rdi
              vmovapd        (%rax), er0
              vmovapd        (%rbx), er2
              vmovapd        (%rcx), er1
              vmovapd        (%rdx), er3
              vpermpd        prmt_cntrl, (%rax, %rdi), ei0
              vpermpd        prmt_cntrl, (%rbx, %rdi), ei2
              vpermpd        prmt_cntrl, (%rcx, %rdi), ei1
              vpermpd        prmt_cntrl, (%rdx, %rdi), ei3
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
              vpermpd        prmt_cntrl, er1, er1
              vpermpd        prmt_cntrl, er2, er2
              vpermpd        prmt_cntrl, ei3, ei3
              vpermpd        prmt_cntrl, ei0, ei0
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
done0:        pop            %r15
              pop            %r14
              pop            %r13
              pop            %r12
              pop            %rbx
              ret



























