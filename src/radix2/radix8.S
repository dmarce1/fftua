#include      "common.h"


              .global        radix8


              .text


radix8:       mov            N, NLO
              shr            $3, NLO
              cmp            $0, k2
              je             keq0
              mov            N2, %rax
              shr            %rax
              cmp            k2, %rax
              je             keqN2o2
              jmp            k2rest
keq0:         hlt
keqN2o2:      hlt
k2rest:       mov            N2, %rdx
              dec            %rdx
              mov            X, %rdi
              sub            X0, %rdi
              shr            $6, %rdi
              bsf            NLO, %rcx
              shr            %rcx, %rdi
              push           %rdi
              mov            N2, %rsi
              call           reverse_neg
              mov            %rax, %rdi
              pop            %rax
              sub            %rax, %rdi
              jl             rcalls
              imul           NLO, %rdi
              shl            $3, %rdi
              mov            N2, %rax
              shr            %rax
              push           X
              push           k2
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
              vmovq          (C, %rax, 8), %xmm0
              vmovq          (S, %rax, 8), %xmm1
              vmulsd         %xmm1, %xmm1, %xmm2
              vmulsd         %xmm0, %xmm1, %xmm3
              vfmsub231sd    %xmm0, %xmm0, %xmm2
              vfmadd231sd    %xmm1, %xmm0, %xmm3
              vmulsd         %xmm1, %xmm3, %xmm4
              vmulsd         %xmm0, %xmm3, %xmm5
              vfmsub231sd    %xmm0, %xmm2, %xmm4
              vfmadd231sd    %xmm1, %xmm2, %xmm5
              vmulsd         %xmm3, %xmm3, %xmm6
              vmulsd         %xmm2, %xmm3, %xmm7
              vfmsub231sd    %xmm2, %xmm2, %xmm6
              vfmadd231sd    %xmm3, %xmm2, %xmm7
              vmulsd         %xmm3, %xmm5, %xmm8
              vmulsd         %xmm2, %xmm5, %xmm9
              vfmsub231sd    %xmm2, %xmm4, %xmm8
              vfmadd231sd    %xmm3, %xmm4, %xmm9
              vmulsd         %xmm5, %xmm5, %xmm10
              vmulsd         %xmm4, %xmm5, %xmm11
              vfmsub231sd    %xmm4, %xmm4, %xmm10
              vfmadd231sd    %xmm5, %xmm4, %xmm11
              vmulsd         %xmm5, %xmm7, %xmm12
              vmulsd         %xmm4, %xmm7, %xmm13
              vfmsub231sd    %xmm4, %xmm6, %xmm12
              vfmadd231sd    %xmm5, %xmm6, %xmm13
              vmovq          %xmm0, COS1
              vmovq          %xmm1, SIN1
              vmovq          %xmm2, COS2
              vmovq          %xmm3, SIN2
              vmovq          %xmm4, COS3
              vmovq          %xmm5, SIN3
              vmovq          %xmm6, COS4
              vmovq          %xmm7, SIN4
              vmovq          %xmm8, COS5
              vmovq          %xmm9, SIN5
              vmovq          %xmm10, COS6
              vmovq          %xmm11, SIN6
              vmovq          %xmm12, COS7
              vmovq          %xmm13, SIN7
              xor            ilo, ilo
k2loop:       mov            NLO, %rsi
              shl            %rsi
              lea            (X, ilo, 8), %rax
              lea            (%rax, %rsi, 8), %rbx
              lea            (%rbx, %rsi, 8), %rcx
              lea            (%rcx, %rsi, 8), %rdx
              vmovapd        (%rax), er0
              vmovapd        (%rbx), er1
              vmovapd        (%rcx), er2
              vmovapd        (%rdx), er3
              vmovapd        (%rax, %rdi, 8), ei0
              vmovapd        (%rbx, %rdi, 8), ei1
              vmovapd        (%rcx, %rdi, 8), ei2
              vmovapd        (%rdx, %rdi, 8), ei3
              vbroadcastsd   COS2, tr0
              vbroadcastsd   SIN2, ti0
              vmulpd         ti0, ei1, tr1
              vmulpd         tr0, ei1, ti1
              vfmsub231pd    tr0, er1, tr1
              vfmadd231pd    ti0, er1, ti1
              vmovapd        ti1, ei1
              vmovapd        tr1, er1
              vbroadcastsd   COS4, tr0
              vbroadcastsd   SIN4, ti0
              vmulpd         ti0, ei2, tr1
              vmulpd         tr0, ei2, ti1
              vfmsub231pd    tr0, er2, tr1
              vfmadd231pd    ti0, er2, ti1
              vmovapd        ti1, ei2
              vmovapd        tr1, er2
              vbroadcastsd   COS6, tr0
              vbroadcastsd   SIN6, ti0
              vmulpd         ti0, ei3, tr1
              vmulpd         tr0, ei3, ti1
              vfmsub231pd    tr0, er3, tr1
              vfmadd231pd    ti0, er3, ti1
              vmovapd        ti1, ei3
              vmovapd        tr1, er3
              vaddpd         er2, er0, tr0
              vaddpd         ei2, ei0, ti0
              vsubpd         er2, er0, tr2
              vsubpd         ei2, ei0, ti2
              vaddpd         er3, er1, tr1
              vaddpd         ei3, ei1, ti1
              vsubpd         er3, er1, tr3
              vsubpd         ei3, ei1, ti3
              vaddpd         tr1, tr0, er0
              vaddpd         ti1, ti0, ei0
              vaddpd         ti3, tr2, er1
              vsubpd         tr3, ti2, ei1
              vsubpd         tr1, tr0, er2
              vsubpd         ti1, ti0, ei2
              vsubpd         ti3, tr2, er3
              vaddpd         tr3, ti2, ei3
              vmovupd        er0, UR4
              vmovupd        er1, UR5
              vmovupd        er2, UR6
              vmovupd        er3, UR7
              vmovupd        ei0, UI4
              vmovupd        ei1, UI5
              vmovupd        ei2, UI6
              vmovupd        ei3, UI7
              lea            (%rax, NLO, 8), %rax
              lea            (%rbx, NLO, 8), %rbx
              lea            (%rcx, NLO, 8), %rcx
              lea            (%rdx, NLO, 8), %rdx
              vmovapd        (%rax), er0
              vmovapd        (%rbx), er1
              vmovapd        (%rcx), er2
              vmovapd        (%rdx), er3
              vmovapd        (%rax, %rdi, 8), ei0
              vmovapd        (%rbx, %rdi, 8), ei1
              vmovapd        (%rcx, %rdi, 8), ei2
              vmovapd        (%rdx, %rdi, 8), ei3
              vbroadcastsd   COS1, tr0
              vbroadcastsd   SIN1, ti0
              vmulpd         ti0, ei0, tr1
              vmulpd         tr0, ei0, ti1
              vfmsub231pd    tr0, er0, tr1
              vfmadd231pd    ti0, er0, ti1
              vmovapd        ti1, ei0
              vmovapd        tr1, er0
              vbroadcastsd   COS3, tr0
              vbroadcastsd   SIN3, ti0
              vmulpd         ti0, ei1, tr1
              vmulpd         tr0, ei1, ti1
              vfmsub231pd    tr0, er1, tr1
              vfmadd231pd    ti0, er1, ti1
              vmovapd        ti1, ei1
              vmovapd        tr1, er1
              vbroadcastsd   COS5, tr0
              vbroadcastsd   SIN5, ti0
              vmulpd         ti0, ei2, tr1
              vmulpd         tr0, ei2, ti1
              vfmsub231pd    tr0, er2, tr1
              vfmadd231pd    ti0, er2, ti1
              vmovapd        ti1, ei2
              vmovapd        tr1, er2
              vbroadcastsd   COS7, tr0
              vbroadcastsd   SIN7, ti0
              vmulpd         ti0, ei3, tr1
              vmulpd         tr0, ei3, ti1
              vfmsub231pd    tr0, er3, tr1
              vfmadd231pd    ti0, er3, ti1
              vmovapd        ti1, ei3
              vmovapd        tr1, er3
              vaddpd         er2, er0, tr0
              vaddpd         ei2, ei0, ti0
              vsubpd         er2, er0, tr2
              vsubpd         ei2, ei0, ti2
              vaddpd         er3, er1, tr1
              vaddpd         ei3, ei1, ti1
              vsubpd         er3, er1, tr3
              vsubpd         ei3, ei1, ti3
              vaddpd         tr1, tr0, er0
              vaddpd         ti1, ti0, ei0
              vaddpd         ti3, tr2, er1
              vsubpd         tr3, ti2, ei1
              vsubpd         tr1, tr0, er2
              vsubpd         ti1, ti0, ei2
              vsubpd         ti3, tr2, er3
              vaddpd         tr3, ti2, ei3
              vmovapd        ti0, ei0
              vmovapd        tr0, er0
              vmulpd         sin45, ei1, tr1
              vmulpd         cos45, ei1, ti1
              vfmsub231pd    cos45, er1, tr1
              vfmadd231pd    sin45, er1, ti1
              vmovapd        ti1, ei1
              vmovapd        tr1, er1
              vmulpd         sin135, ei3, tr1
              vmulpd         cos135, ei3, ti1
              vfmsub231pd    cos135, er3, tr1
              vfmadd231pd    sin135, er3, ti1
              vmovapd        ti1, ei3
              vmovapd        tr1, er3
              vmovapd        ei2, tr2
              vmulpd         none, er2, ti2
              vmovapd        tr2, er2
              vmovapd        ti2, ei2
              vaddpd         UR4, er0, tr0
              vaddpd         UI4, ei0, ti0
              vsubpd         UR4, er0, er0
              vsubpd         UI4, ei0, ei0
              vaddpd         UR5, er1, tr1
              vaddpd         UI5, ei1, ti1
              vsubpd         UR5, er1, er1
              vsubpd         UI5, ei1, ei1
              vaddpd         UR6, er2, tr2
              vaddpd         UI6, ei2, ti2
              vsubpd         UR6, er2, er2
              vsubpd         UI6, ei2, ei2
              vaddpd         UR7, er3, tr3
              vaddpd         UI7, ei3, ti3
              vsubpd         UR7, er3, er3
              vsubpd         UI7, ei3, ei3
              vmulpd         none, er0, er0
              vmulpd         none, er1, er1
              vmulpd         none, er2, er2
              vmulpd         none, er3, er3
              vmovapd        ei3, (%rax)
              vmovapd        ei2, (%rbx)
              vmovapd        ei1, (%rcx)
              vmovapd        ei0, (%rdx)
              vmovapd        ti3, (%rax, %rdi, 8)
              vmovapd        ti2, (%rbx, %rdi, 8)
              vmovapd        ti1, (%rcx, %rdi, 8)
              vmovapd        ti0, (%rdx, %rdi, 8)
              mov            NLO, %rsi
              neg            %rsi
              lea            (%rax, %rsi, 8), %rax
              lea            (%rbx, %rsi, 8), %rbx
              lea            (%rcx, %rsi, 8), %rcx
              lea            (%rdx, %rsi, 8), %rdx
              vmovapd        tr0, (%rax)
              vmovapd        tr1, (%rbx)
              vmovapd        tr2, (%rcx)
              vmovapd        tr3, (%rdx)
              vmovapd        er0, (%rax, %rdi, 8)
              vmovapd        er1, (%rbx, %rdi, 8)
              vmovapd        er2, (%rcx, %rdi, 8)
              vmovapd        er3, (%rdx, %rdi, 8)
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
              mov            NLO, %rdx
              shl            %rdx
              neg            %rdx
              mov            $7, %rax
              mov            $7, %rbx
              shl            %rcx, %rax
              or             k2, %rax
              lea            (X, %rbx, 8), %rcx
              push           %rcx
              push           %rax
              bsf            N2, %rcx
              mov            NLO, %rdx
              shl            %rdx
              neg            %rdx
              mov            $6, %rax
              mov            $3, %rbx
              shl            %rcx, %rax
              or             k2, %rax
              lea            (X, %rbx, 8), %rcx
              push           %rcx
              push           %rax
              bsf            N2, %rcx
              mov            NLO, %rdx
              shl            %rdx
              neg            %rdx
              mov            $5, %rax
              mov            $5, %rbx
              shl            %rcx, %rax
              or             k2, %rax
              lea            (X, %rbx, 8), %rcx
              push           %rcx
              push           %rax
              bsf            N2, %rcx
              mov            NLO, %rdx
              shl            %rdx
              neg            %rdx
              mov            $4, %rax
              mov            $1, %rbx
              shl            %rcx, %rax
              or             k2, %rax
              lea            (X, %rbx, 8), %rcx
              push           %rcx
              push           %rax
              bsf            N2, %rcx
              mov            NLO, %rdx
              shl            %rdx
              neg            %rdx
              mov            $3, %rax
              mov            $6, %rbx
              shl            %rcx, %rax
              or             k2, %rax
              lea            (X, %rbx, 8), %rcx
              push           %rcx
              push           %rax
              bsf            N2, %rcx
              mov            NLO, %rdx
              shl            %rdx
              neg            %rdx
              mov            $2, %rax
              mov            $2, %rbx
              shl            %rcx, %rax
              or             k2, %rax
              lea            (X, %rbx, 8), %rcx
              push           %rcx
              push           %rax
              bsf            N2, %rcx
              mov            NLO, %rdx
              shl            %rdx
              neg            %rdx
              mov            $1, %rax
              mov            $4, %rbx
              shl            %rcx, %rax
              or             k2, %rax
              lea            (X, %rbx, 8), %rcx
              push           %rcx
              push           %rax
              bsf            N2, %rcx
              mov            NLO, %rdx
              shl            %rdx
              neg            %rdx
              mov            $0, %rax
              mov            $0, %rbx
              shl            %rcx, %rax
              or             k2, %rax
              lea            (X, %rbx, 8), %rcx
              mov            %rcx, X
              mov            %rax, k2
              mov            N, %rax
              shr            $3, %rax
              mov            %rax, N
              shl            $3, N2
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
              shl            $3, %rax
              mov            %rax, N
              shr            $3, N2
done1:        ret

