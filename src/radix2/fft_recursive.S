#include      "common.h"


              .global        fft_recursive


              .text

fft_recursive:enter          STACK_SIZE, $0
              push           %rbx
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
              movsd          tcos2, (%rbx)
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
              movsd          ur2, tcos2
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
              movsd          8(C), tcos1
              movsd          8(S), tsin1
              vmulsd         tsin1, tsin1, tcos2
              vmulsd         tcos1, tsin1, tsin2
              vfmsub231sd    tcos1, tcos1, tcos2
              vfmadd231sd    tsin1, tcos1, tsin2
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
              vmovapd        two, ytwo
              vfmadd231sd    tsin2, ui2, sr0
              vfnmadd231sd   tsin2, ur2, si0
              vfnmadd132sd   tcos2, sr0, ur2
              vfnmadd132sd   tcos2, si0, ui2
              vfmsub132sd    xtwo, ur2, ur0
              vfmsub132sd    xtwo, ui2, ui0
              movsd          ur1, sr0
              movsd          ui1, si0
              vfmadd231sd    tsin2, ui3, sr0
              vfnmadd231sd   tsin2, ur3, si0
              vfnmadd132sd   tcos2, sr0, ur3
              vfnmadd132sd   tcos2, si0, ui3
              vfmsub132sd    xtwo, ur3, ur1
              vfmsub132sd    xtwo, ui3, ui1
              movsd          ur0, sr0
              movsd          ui0, si0
              vfmadd231sd    tsin1, ui1, sr0
              vfnmadd231sd   tsin1, ur1, si0
              vfnmadd132sd   tcos1, sr0, ur1
              vfmsub132sd    tcos1, si0, ui1
              vfmsub132sd    xtwo, ur1, ur0
              vfmadd132sd    xtwo, ui1, ui0
              movsd          ur2, sr0
              movsd          ui2, si0
              vfmadd231sd    tcos1, ui3, sr0
              vfnmadd231sd   tcos1, ur3, si0
              vfmadd132sd    tsin1, sr0, ur3
              vfmadd132sd    tsin1, si0, ui3
              vfmsub132sd    xtwo, ur3, ur2
              vfnmadd132sd   xtwo, ui3, ui2
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
              movlpd         16(C), tcos1
              movlpd         16(S), tsin1
              vmovhpd        24(C), tcos1, tcos1
              vmovhpd        24(S), tsin1, tsin1
              vmulpd         sin1, sin1, cos2
              vmulpd         cos1, sin1, sin2
              vfmsub231pd    cos1, cos1, cos2
              vfmadd231pd    sin1, cos1, sin2
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
              vfmadd231pd    tsin2, ui2, sr0
              vfnmadd231pd   tsin2, ur2, si0
              vfnmadd132pd   tcos2, sr0, ur2
              vfnmadd132pd   tcos2, si0, ui2
              vfmsub132pd    xtwo, ur2, ur0
              vfmsub132pd    xtwo, ui2, ui0
              vmovapd        ur1, sr0
              vmovapd        ui1, si0
              vfmadd231pd    tsin2, ui3, sr0
              vfnmadd231pd   tsin2, ur3, si0
              vfnmadd132pd   tcos2, sr0, ur3
              vfnmadd132pd   tcos2, si0, ui3
              vfmsub132pd    xtwo, ur3, ur1
              vfmsub132pd    xtwo, ui3, ui1
              vmovapd        ur0, sr0
              vmovapd        ui0, si0
              vfmadd231pd    tsin1, ui1, sr0
              vfnmadd231pd   tsin1, ur1, si0
              vfnmadd132pd   tcos1, sr0, ur1
              vfmsub132pd    tcos1, si0, ui1
              vfmsub132pd    xtwo, ur1, ur0
              vfmadd132pd    xtwo, ui1, ui0
              vmovapd        ur2, sr0
              vmovapd        ui2, si0
              vfmadd231pd    tcos1, ui3, sr0
              vfnmadd231pd   tcos1, ur3, si0
              vfmadd132pd    tsin1, sr0, ur3
              vfmadd132pd    tsin1, si0, ui3
              vfmsub132pd    xtwo, ur3, ur2
              vfnmadd132pd   xtwo, ui3, ui2
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
final:        vmovapd        (C, k2, 8), cos1
              vmovapd        (S, k2, 8), sin1
              vmulpd         sin1, sin1, cos2
              vmulpd         cos1, sin1, sin2
              vfmsub231pd    cos1, cos1, cos2
              vfmadd231pd    sin1, cos1, sin2
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
              vfmadd231pd    sin2, ei2, tr0
              vfnmadd231pd   sin2, er2, ti0
              vfnmadd132pd   cos2, tr0, er2
              vfnmadd132pd   cos2, ti0, ei2
              vfmsub132pd    ytwo, er2, er0
              vfmsub132pd    ytwo, ei2, ei0
              vmovapd        er1, tr0
              vmovapd        ei1, ti0
              vfmadd231pd    sin2, ei3, tr0
              vfnmadd231pd   sin2, er3, ti0
              vfnmadd132pd   cos2, tr0, er3
              vfnmadd132pd   cos2, ti0, ei3
              vfmsub132pd    ytwo, er3, er1
              vfmsub132pd    ytwo, ei3, ei1
              vmovapd        er0, tr0
              vmovapd        ei0, ti0
              vfmadd231pd    sin1, ei1, tr0
              vfnmadd231pd   sin1, er1, ti0
              vfnmadd132pd   cos1, tr0, er1
              vfmsub132pd    cos1, ti0, ei1
              vfmsub132pd    ytwo, er1, er0
              vfmadd132pd    ytwo, ei1, ei0
              vmovapd        er2, tr0
              vmovapd        ei2, ti0
              vfmadd231pd    cos1, ei3, tr0
              vfnmadd231pd   cos1, er3, ti0
              vfmadd132pd    sin1, tr0, er3
              vfmadd132pd    sin1, ti0, ei3
              vfmsub132pd    ytwo, er3, er2
              vfnmadd132pd   ytwo, ei3, ei2
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
              leave
              ret



























