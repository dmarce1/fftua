#include      "common.h"

              .global        radix2

              .text


radix2:       mov            N, NLO
              shr            NLO
              cmp            $0, k2
              je             r2keq0
              mov            N2, %rax
              shr            %rax
              cmp            k2, %rax
              je             r2keqN2o2
              jmp            r2k2rest
r2keq0:       xor            ilo, ilo
r2keq0lp:     lea            (X, ilo, 8), %rax
              lea            (%rax, NLO, 8), %rbx
              vmovapd        (%rax), er0
              vmovapd        (%rbx), er2
              vaddpd         er2, er0, tr0
              vsubpd         er2, er0, tr2
              vmovapd        tr0, (%rax)
              vmovapd        tr2, (%rbx)
              add            $4, ilo
              cmp            ilo, NLO
              jg             r2keq0lp
              jmp            r2rcalls
r2keqN2o2:    xor            ilo, ilo
r2N2o2lp:     lea            (X, ilo, 8), %rax
              lea            (%rax, NLO, 8), %rbx
              vmovapd        (%rbx), er1
              vmulpd         none, er1, er1
              vmovapd        er1, (%rbx)
              add            $4, ilo
              cmp            ilo, NLO
              jg             r2N2o2lp
              jmp            r2rcalls
r2k2rest:     mov            X, %rdi
              sub            X0, %rdi
              shr            $3, %rdi
              bsf            NLO, %rcx
              shr            %rcx, %rdi
              mov            %rdi, %rsi
              bsr            %rdi, %rcx
              mov            $3, %rdi
              shl            %rcx, %rdi
              dec            %rdi
              sub            %rsi, %rdi
              sub            %rsi, %rdi
              jl             r2rcalls
              imul           NLO, %rdi
              shl            %rdi
              mov            N2, %rax
              shr            %rax
              push           X
              push           k2
              cmp            k2, %rax
              jg             r2skipneg
              lea            (X, %rdi, 8), X
              neg            %rdi
              mov            N2, %rax
              sub            k2, %rax
              mov            %rax, k2
r2skipneg:    mov            k2, %rax
              imul           NLO, %rax
              vbroadcastsd   (C, %rax, 8), tr1
              vbroadcastsd   (S, %rax, 8), ti1
              xor            ilo, ilo
r2k2loop:     lea            (X, ilo, 8), %rax
              lea            (%rax, NLO, 8), %rbx
              vmovapd        (%rax), er0
              vmovapd        (%rbx), er1
              vmovapd        (%rax, %rdi, 8), ei0
              vmovapd        (%rbx, %rdi, 8), ei1
              vmovapd        er0, tr0
              vmovapd        ei0, ti0
              vfmadd231pd    ti1, ei1, tr0
              vfnmadd231pd   ti1, er1, ti0
              vfnmadd132pd   tr1, tr0, er1
              vfmsub132pd    tr1, ti0, ei1
              vfmsub132pd    two, er1, er0
              vfmadd132pd    two, ei1, ei0
              vmovapd        er0, (%rax)
              vmovapd        ei1, (%rbx)
              vmovapd        er1, (%rax, %rdi, 8)
              vmovapd        ei0, (%rbx, %rdi, 8)
              add            $4, ilo
              cmp            ilo, NLO
              jg             r2k2loop
              pop            k2
              pop            X
r2rcalls:     cmp            $4, NLO
              jle            r2done1
              push           k2
              push           X
              bsf            N2, %rcx
              mov            $1, %rax
              shl            %rcx, %rax
              or             k2, %rax
              push           %rax
              lea            (X, NLO, 8), %rax
              push           %rax
              mov            N, %rax
              shr            %rax
              mov            %rax, N
              shl            N2
              call           fft_dispatch
              pop            X
              pop            k2
              call           fft_dispatch
              pop            X
              pop            k2
              mov            N, %rax
              shl            %rax
              mov            %rax, N
              shr            N2
              mov            N, NLO
              shr            NLO
r2done1:      ret


