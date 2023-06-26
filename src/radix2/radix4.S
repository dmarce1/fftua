#include      "common.h"


              .global        radix4


              .text


radix4:       mov            N, NLO
              shr            $2, NLO
              cmp            $0, k2
              je             keq0
              mov            N2, %rax
              shr            %rax
              cmp            k2, %rax
              je             keqN2o2
              jmp            k2rest
keq0:         xor            ilo, ilo
keq0lp:       lea            (X, ilo, 8), %rax
              lea            (%rax, NLO, 8), %rbx
              lea            (%rbx, NLO, 8), %rcx
              lea            (%rcx, NLO, 8), %rdx
              vmovapd        (%rax), er0
              vmovapd        (%rbx), er1
              vmovapd        (%rcx), er2
              vmovapd        (%rdx), er3
              vaddpd         er2, er0, tr0
              vsubpd         er2, er0, tr2
              vaddpd         er3, er1, tr1
              vsubpd         er1, er3, tr3
              vaddpd         tr1, tr0, er0
              vsubpd         tr1, tr0, er2
              vmovapd        er0, (%rax)
              vmovapd        er2, (%rbx)
              vmovapd        tr2, (%rcx)
              vmovapd        tr3, (%rdx)
              add            $4, ilo
              cmp            ilo, NLO
              jg             keq0lp
              jmp            rcalls
keqN2o2:      xor            ilo, ilo
N2o2lp:       lea            (X, ilo, 8), %rax
              lea            (%rax, NLO, 8), %rbx
              lea            (%rbx, NLO, 8), %rcx
              lea            (%rcx, NLO, 8), %rdx
              vmovapd        (%rax), er0
              vmovapd        (%rbx), er1
              vmovapd        (%rcx), er2
              vmovapd        (%rdx), er3
              vaddpd         er3, er1, tr0
              vsubpd         er3, er1, tr2
              vmulpd         tw45, tr2, tr1
              vmulpd         tw45, tr0, tr3
              vmovapd        er0, tr0
              vmovapd        er2, tr2
              vaddpd         tr1, tr0, er0
              vaddpd         tr3, tr2, er3
              vsubpd         tr1, tr0, er1
              vsubpd         tr3, tr2, er2
              vmulpd         none, er3, er3
              vmovapd        er0, (%rax)
              vmovapd        er2, (%rbx)
              vmovapd        er1, (%rcx)
              vmovapd        er3, (%rdx)
              add            $4, ilo
              cmp            ilo, NLO
              jg             N2o2lp
              jmp            rcalls
k2rest:       mov            X, %rdi
              sub            X0, %rdi
              shr            $5, %rdi
              bsf            NLO, %rcx
              shr            %rcx, %rdi
              mov            %rdi, %rsi
              bsr            %rdi, %rcx
              mov            $3, %rdi
              shl            %rcx, %rdi
              dec            %rdi
              sub            %rsi, %rdi
              sub            %rsi, %rdi
              jl             rcalls
              imul           NLO, %rdi
              shl            $2, %rdi
              push           X
              push           k2
              mov            N2, %rax
              shr            %rax
              cmp            k2, %rax
              jg             skipneg
              lea            (X, %rdi, 8), X
              neg            %rdi
              mov            N2, %rax
              sub            k2, %rax
              mov            %rax, k2
skipneg:      vmovapd        two, ytwo
		      mov            k2, %rax
              imul           NLO, %rax
              vbroadcastsd   (C, %rax, 8), cos1
              vbroadcastsd   (S, %rax, 8), sin1
              vmulpd         sin1, sin1, cos2
              vmulpd         cos1, sin1, sin2
              vfmsub231pd    cos1, cos1, cos2
              vfmadd231pd    sin1, cos1, sin2
              xor            ilo, ilo
k2loop:       lea            (X, ilo, 8), %rax
              lea            (%rax, NLO, 8), %rbx
              lea            (%rbx, NLO, 8), %rcx
              lea            (%rcx, NLO, 8), %rdx
              vmovapd        (%rax), er0
              vmovapd        (%rbx), er1
              vmovapd        (%rcx), er2
              vmovapd        (%rdx), er3
              vmovapd        (%rax, %rdi, 8), ei0
              vmovapd        (%rbx, %rdi, 8), ei1
              vmovapd        (%rcx, %rdi, 8), ei2
              vmovapd        (%rdx, %rdi, 8), ei3
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
              vmovapd        er0, (%rax)
              vmovapd        ei1, (%rbx)
              vmovapd        er3, (%rcx)
              vmovapd        ei2, (%rdx)
              vmovapd        er2, (%rax, %rdi, 8)
              vmovapd        ei3, (%rbx, %rdi, 8)
              vmovapd        er1, (%rcx, %rdi, 8)
              vmovapd        ei0, (%rdx, %rdi, 8)
              add            $4, ilo
              cmp            ilo, NLO
              jg             k2loop
              pop            k2
              pop            X
rcalls:       cmp            $4, NLO
              jle            done1
              push           X
              push           k2
              bsf            N2, %rcx
              imul           $3, NLO, %rbx
              mov            NLO, %rdx
              shl            %rdx
              neg            %rdx
              mov            $3, %rax
              shl            %rcx, %rax
              or             k2, %rax
              lea            (X, %rbx, 8), X
              push           X
              push           %rax
              mov            $2, %rax
              shl            %rcx, %rax
              or             k2, %rax
              lea            (X, %rdx, 8), X
              push           X
              push           %rax
              mov            $1, %rax
              shl            %rcx, %rax
              or             k2, %rax
              lea            (X, NLO, 8), X
              push           X
              push           %rax
              lea            (X, %rdx, 8), X
              mov            N, %rax
              shr            $2, %rax
              mov            %rax, N
              shl            $2, N2
              call           fft_dispatch
              pop            k2
              pop            X
              call           fft_dispatch
              pop            k2
              pop            X
              call           fft_dispatch
              pop            k2
              pop            X
              call           fft_dispatch
              pop            k2
              pop            X
              mov            N, %rax
              shl            $2, %rax
              mov            %rax, N
              shr            $2, N2
done1:        ret

