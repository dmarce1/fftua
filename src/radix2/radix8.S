#include      "common.h"


              .global        radix8


              .align         32
C0:           .double        3.82683432365089782e-01
              .double        3.82683432365089782e-01
              .double        3.82683432365089782e-01
              .double        3.82683432365089782e-01
C1:           .double        9.23879532511286738e-01
              .double        9.23879532511286738e-01
              .double        9.23879532511286738e-01
              .double        9.23879532511286738e-01
C2:           .double        7.07106781186547573e-01
              .double        7.07106781186547573e-01
              .double        7.07106781186547573e-01
              .double        7.07106781186547573e-01


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
keq0:         xor            ilo, ilo
              mov            NLO, %rdi
              shl            $3, %rdi
keq0loop:     lea            (X, ilo, 8), %rax
              lea            (%rax, %rdi, 2), %rbx
              lea            (%rbx, %rdi, 2), %rcx
              lea            (%rcx, %rdi, 2), %rdx
              vmovapd        (%rax), %ymm0
              vmovapd        (%rax, %rdi, 1), %ymm1
              vmovapd        (%rbx), %ymm2
              vmovapd        (%rbx, %rdi, 1), %ymm3
              vmovapd        (%rcx), %ymm4
              vmovapd        (%rcx, %rdi, 1), %ymm5
              vmovapd        (%rdx), %ymm6
              vmovapd        (%rdx, %rdi, 1), %ymm7
	          vaddpd         %ymm0, %ymm4, %ymm8
	          vsubpd         %ymm4, %ymm0, %ymm0
	          vaddpd         %ymm2, %ymm6, %ymm4
	          vsubpd         %ymm6, %ymm2, %ymm2
	          vaddpd         %ymm8, %ymm4, %ymm6
	          vmovapd        %ymm2, %ymm9
	          vsubpd         %ymm4, %ymm8, %ymm2
	          vaddpd         %ymm1, %ymm5, %ymm8
	          vsubpd         %ymm5, %ymm1, %ymm1
	          vaddpd         %ymm7, %ymm3, %ymm4
	          vsubpd         %ymm3, %ymm7, %ymm7
	          vsubpd         %ymm4, %ymm8, %ymm3
	          vaddpd         %ymm8, %ymm4, %ymm8
	          vmovapd        %ymm0, %ymm4
	          vaddpd         %ymm6, %ymm8, %ymm0
	          vmovapd        %ymm6, %ymm5
	          vmulpd         none, %ymm3, %ymm6
	          vmovapd        %ymm4, %ymm3
	          vsubpd         %ymm8, %ymm5, %ymm4
	          vmulpd         tw45, %ymm7, %ymm7
	          vmulpd         tw45, %ymm1, %ymm1
	          vaddpd         %ymm1, %ymm7, %ymm8
	          vmovapd        %ymm1, %ymm5
	          vaddpd         %ymm3, %ymm8,%ymm1
	          vsubpd         %ymm5, %ymm7, %ymm7
	          vmovapd        %ymm7, %ymm5
	          vsubpd         %ymm9, %ymm5, %ymm7
	          vmovapd        %ymm3, %ymm10
	          vsubpd         %ymm8, %ymm10, %ymm3
	          vmovapd        %ymm5, %ymm8
	          vaddpd         %ymm8, %ymm9, %ymm5
              vmovapd        %ymm0, (%rax)
              vmovapd        %ymm4, (%rax, %rdi, 1)
              vmovapd        %ymm2, (%rbx)
              vmovapd        %ymm6, (%rbx, %rdi, 1)
              vmovapd        %ymm1, (%rcx)
              vmovapd        %ymm5, (%rcx, %rdi, 1)
              vmovapd        %ymm3, (%rdx)
              vmovapd        %ymm7, (%rdx, %rdi, 1)
              add            $4, ilo
              cmp            ilo, NLO
              jg             keq0loop
		      jmp            rcalls
keqN2o2:      xor            ilo, ilo
              mov            NLO, %rdi
              shl            $3, %rdi
keqN2o2loop:  lea            (X, ilo, 8), %rax
              lea            (%rax, %rdi, 2), %rbx
              lea            (%rbx, %rdi, 2), %rcx
              lea            (%rcx, %rdi, 2), %rdx
              vmovapd        (%rax), %ymm0
              vmovapd        (%rax, %rdi, 1), %ymm1
              vmovapd        (%rbx), %ymm2
              vmovapd        (%rbx, %rdi, 1), %ymm3
              vmovapd        (%rcx), %ymm4
              vmovapd        (%rcx, %rdi, 1), %ymm5
              vmovapd        (%rdx), %ymm6
              vmovapd        (%rdx, %rdi, 1), %ymm7
	          vmulpd         C0, %ymm5, %ymm8
	          vmulpd         C1, %ymm1, %ymm9
	          vsubpd         %ymm8, %ymm9, %ymm9
	          vmulpd         C1, %ymm5, %ymm5
	          vmulpd         C0, %ymm1, %ymm1
	          vaddpd         %ymm5, %ymm1, %ymm1
	          vmulpd         C2, %ymm6, %ymm6
	          vmulpd         C2, %ymm2, %ymm2
	          vsubpd         %ymm6, %ymm2, %ymm8
	          vaddpd         %ymm6, %ymm2, %ymm2
	          vmulpd         C1, %ymm7, %ymm5
	          vmulpd         C0, %ymm3, %ymm6
	          vsubpd         %ymm5, %ymm6, %ymm6
	          vmulpd         C0, %ymm7, %ymm7
	          vmulpd         C1, %ymm3, %ymm3
	          vaddpd         %ymm7, %ymm3, %ymm3
	          vaddpd         %ymm0, %ymm8, %ymm5
	          vaddpd         %ymm4, %ymm2, %ymm7
	          vsubpd         %ymm8, %ymm0, %ymm0
	          vsubpd         %ymm4, %ymm2, %ymm2
	          vaddpd         %ymm9, %ymm6, %ymm8
	          vaddpd         %ymm1, %ymm3, %ymm4
	          vsubpd         %ymm6, %ymm9, %ymm9
	          vsubpd         %ymm1, %ymm3, %ymm3
	          vmovapd        %ymm0, %ymm1
	          vaddpd         %ymm5, %ymm8, %ymm0
	          vaddpd         %ymm7, %ymm4, %ymm6
	          vmovapd        %ymm7,%ymm10
	          vmulpd         none, %ymm6, %ymm7
	          vmovapd        %ymm3, %ymm6
	          vsubpd         %ymm8, %ymm5, %ymm3
	          vsubpd         %ymm10, %ymm4, %ymm4
	          vmovapd        %ymm2, %ymm8
	          vaddpd         %ymm1, %ymm6, %ymm2
	          vsubpd         %ymm9, %ymm8, %ymm5
	          vmovapd        %ymm1, %ymm10
	          vsubpd         %ymm6, %ymm10, %ymm1
	          vaddpd         %ymm8, %ymm9, %ymm9
	          vmulpd         none, %ymm9, %ymm6
	          vmovapd        %ymm4, %ymm8
	          vmulpd         none, %ymm8, %ymm4
              vmovapd        %ymm0, (%rax)
              vmovapd        %ymm4, (%rax, %rdi, 1)
              vmovapd        %ymm2, (%rbx)
              vmovapd        %ymm6, (%rbx, %rdi, 1)
              vmovapd        %ymm1, (%rcx)
              vmovapd        %ymm5, (%rcx, %rdi, 1)
              vmovapd        %ymm3, (%rdx)
              vmovapd        %ymm7, (%rdx, %rdi, 1)
              add            $4, ilo
              cmp            ilo, NLO
              jg             keqN2o2loop
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
skipneg:      mov            k2, %rax
              imul           NLO, %rax
              imul           $2, %rax, %rbx
              imul           $3, %rax, %rcx
              imul           $4, %rax, %rdx
              vbroadcastsd   (C, %rax, 8), %ymm2
              vbroadcastsd   (S, %rax, 8), %ymm3
              vbroadcastsd   (C, %rbx, 8), %ymm4
              vbroadcastsd   (S, %rbx, 8), %ymm5
              vbroadcastsd   (C, %rcx, 8), %ymm6
              vbroadcastsd   (S, %rcx, 8), %ymm7
              vbroadcastsd   (C, %rdx, 8), %ymm8
              vbroadcastsd   (S, %rdx, 8), %ymm9
              imul           $5, %rax, %rbx
              imul           $6, %rax, %rcx
              imul           $7, %rax, %rdx
              vbroadcastsd   (C, %rbx, 8), %ymm10
              vbroadcastsd   (S, %rbx, 8), %ymm11
              vbroadcastsd   (C, %rcx, 8), %ymm12
              vbroadcastsd   (S, %rcx, 8), %ymm13
              vbroadcastsd   (C, %rdx, 8), %ymm14
              vbroadcastsd   (S, %rdx, 8), %ymm15
              vmovupd        %ymm2, COS1
              vmovupd        %ymm3, SIN1
              vmovupd        %ymm4, COS2
              vmovupd        %ymm5, SIN2
              vmovupd        %ymm6, COS3
              vmovupd        %ymm7, SIN3
              vmovupd        %ymm8, COS4
              vmovupd        %ymm9, SIN4
              vmovupd        %ymm10, COS5
              vmovupd        %ymm11, SIN5
              vmovupd        %ymm12, COS6
              vmovupd        %ymm13, SIN6
              vmovupd        %ymm14, COS7
              vmovupd        %ymm15, SIN7
              xor            ilo, ilo
              mov            NLO, %rsi
              shl            $3, %rsi
k2loop:       lea            (X, ilo, 8), %rax
              lea            (%rax, %rsi, 2), %rbx
              lea            (%rbx, %rsi, 2), %rcx
              lea            (%rcx, %rsi, 2), %rdx
              vmovapd        (%rax), %ymm0
              vmovapd        (%rax, %rdi, 1), %ymm1
              vmovapd        (%rbx), %ymm2
              vmovapd        (%rbx, %rdi, 1), %ymm3
              vmovapd        (%rcx), %ymm4
              vmovapd        (%rcx, %rdi, 1), %ymm5
              vmovapd        (%rdx), %ymm6
              vmovapd        (%rdx, %rdi, 1), %ymm7
              push           %rdx
              push           %rcx
              push           %rbx
              push           %rax
              lea            (%rax, %rdi, 8), %rax
              lea            (%rbx, %rdi, 8), %rbx
              lea            (%rcx, %rdi, 8), %rcx
              lea            (%rdx, %rdi, 8), %rdx
              vmovapd        (%rax), %ymm8
              vmovapd        (%rax, %rdi, 1), %ymm9
              vmovapd        (%rbx), %ymm10
              vmovapd        (%rbx, %rdi, 1), %ymm11
              vmovapd        (%rcx), %ymm12
              vmovapd        (%rcx, %rdi, 1), %ymm13
              vmovapd        (%rdx), %ymm14
              vmovapd        (%rdx, %rdi, 1), %ymm15
              vmovupd        %ymm0, UR0
              vmovupd        %ymm1, UI0
              vmulpd         SIN1, %ymm3, %ymm0
              vmulpd         COS1, %ymm3, %ymm1
              vfmsub231pd    COS1, %ymm2, %ymm0
              vfmadd231pd    SIN1, %ymm2, %ymm1
              vmovapd        %ymm0, %ymm2
              vmovapd        %ymm1, %ymm3
              vmulpd         SIN2, %ymm5, %ymm0
              vmulpd         COS2, %ymm5, %ymm1
              vfmsub231pd    COS2, %ymm4, %ymm0
              vfmadd231pd    SIN2, %ymm4, %ymm1
              vmovapd        %ymm0, %ymm4
              vmovapd        %ymm1, %ymm5
              vmulpd         SIN3, %ymm7, %ymm0
              vmulpd         COS3, %ymm7, %ymm1
              vfmsub231pd    COS3, %ymm6, %ymm0
              vfmadd231pd    SIN3, %ymm6, %ymm1
              vmovapd        %ymm0, %ymm6
              vmovapd        %ymm1, %ymm7
              vmulpd         SIN4, %ymm9, %ymm0
              vmulpd         COS4, %ymm9, %ymm1
              vfmsub231pd    COS4, %ymm8, %ymm0
              vfmadd231pd    SIN4, %ymm8, %ymm1
              vmovapd        %ymm0, %ymm8
              vmovapd        %ymm1, %ymm9
              vmulpd         SIN5, %ymm11, %ymm0
              vmulpd         COS5, %ymm11, %ymm1
              vfmsub231pd    COS5, %ymm10, %ymm0
              vfmadd231pd    SIN5, %ymm10, %ymm1
              vmovapd        %ymm0, %ymm10
              vmovapd        %ymm1, %ymm11
              vmulpd         SIN6, %ymm13, %ymm0
              vmulpd         COS6, %ymm13, %ymm1
              vfmsub231pd    COS6, %ymm12, %ymm0
              vfmadd231pd    SIN6, %ymm12, %ymm1
              vmovapd        %ymm0, %ymm12
              vmovapd        %ymm1, %ymm13
              vmulpd         SIN7, %ymm15, %ymm0
              vmulpd         COS7, %ymm15, %ymm1
              vfmsub231pd    COS7, %ymm14, %ymm0
              vfmadd231pd    SIN7, %ymm14, %ymm1
              vmovapd        %ymm0, %ymm14
              vmovapd        %ymm1, %ymm15
              vmovupd        UR0, %ymm0
              vmovupd        UI0, %ymm1
	          vmovupd        %ymm8, UR3
              vmovupd        %ymm9, UI3
	          vaddpd         %ymm0, %ymm8, %ymm8
	          vaddpd         %ymm1, %ymm9, %ymm9
	          vmovupd        %ymm8, UR0
	          vmovupd        %ymm9, UR1
	          vmovupd        UR3, %ymm8
	          vmovupd        UI3, %ymm9
 	          vsubpd         %ymm8, %ymm0, %ymm0
              vsubpd         %ymm9, %ymm1, %ymm1
	          vaddpd         %ymm4, %ymm12, %ymm8
	          vaddpd         %ymm5, %ymm13, %ymm9
	          vsubpd         %ymm12, %ymm4, %ymm4
	          vsubpd         %ymm13, %ymm5, %ymm5
	          vaddpd         UR0, %ymm8, %ymm12
	          vaddpd         UR1, %ymm9, %ymm13
	          vmovupd        %ymm0, UR3
              vmovupd        %ymm1, UI3
              vmovupd        UR0, %ymm0
              vmovupd        UR1, %ymm1
              vsubpd         %ymm8, %ymm0, %ymm0
              vsubpd         %ymm9, %ymm1, %ymm1
              vmovupd        %ymm0, UR0
              vmovupd        %ymm1, UR1
	          vmovupd        UR3, %ymm0
	          vmovupd        UI3, %ymm1
	          vaddpd         %ymm0, %ymm5, %ymm8
	          vsubpd         %ymm4, %ymm1, %ymm9
	          vsubpd         %ymm5, %ymm0, %ymm0
	          vaddpd         %ymm1, %ymm4, %ymm4
	          vaddpd         %ymm2, %ymm10, %ymm1
	          vaddpd         %ymm3, %ymm11, %ymm5
	          vsubpd         %ymm10, %ymm2, %ymm2
	          vsubpd         %ymm11, %ymm3, %ymm3
	          vaddpd         %ymm6, %ymm14, %ymm10
	          vaddpd         %ymm7, %ymm15, %ymm11
	          vsubpd         %ymm14, %ymm6, %ymm6
	          vsubpd         %ymm15, %ymm7, %ymm7
	          vaddpd         %ymm1, %ymm10, %ymm14
	          vaddpd         %ymm5, %ymm11, %ymm15
	          vsubpd         %ymm10, %ymm1, %ymm1
	          vsubpd         %ymm11, %ymm5, %ymm5
	          vaddpd         %ymm7, %ymm2, %ymm10
	          vsubpd         %ymm6, %ymm3, %ymm11
	          vsubpd         %ymm7, %ymm2, %ymm2
	          vaddpd         %ymm3, %ymm6, %ymm6
	          vmulpd         C2, %ymm10, %ymm10
	          vmulpd         C2, %ymm11, %ymm11
	          vaddpd         %ymm10, %ymm11, %ymm3
	          vsubpd         %ymm10, %ymm11, %ymm11
	          vmulpd         C2, %ymm2, %ymm6
	          vmulpd         C2, %ymm2, %ymm6
	          vsubpd         %ymm2, %ymm6, %ymm10
	          vaddpd         %ymm6, %ymm2, %ymm2
	          vmovapd        %ymm0, %ymm6
	          vaddpd         %ymm14, %ymm12, %ymm0
	          vmovapd        %ymm1, %ymm7
	          vaddpd         %ymm15, %ymm13, %ymm1
	          vmovapd        %ymm8, UR2
	          vsubpd         %ymm14, %ymm12, %ymm8
	          vmovapd        %ymm9, %ymm12
	          vsubpd         %ymm15, %ymm13, %ymm9
	          vmovapd        %ymm2, %ymm13
	          vaddpd         UR2, %ymm3, %ymm2
	          vmovapd        %ymm3, %ymm14
	          vaddpd         %ymm11, %ymm12, %ymm3
	          vmovapd        %ymm10, %ymm15
	          vsubpd         UR2, %ymm14, %ymm10
	          vmulpd         none, %ymm10, %ymm10
	          vmovupd        %ymm11, UR2
	          vsubpd         UR2, %ymm12, %ymm11
	          vmovupd        %ymm4, UR2
	          vaddpd         UR0, %ymm5, %ymm4
	          vmovapd        %ymm5, %ymm12
	          vsubpd         UR1, %ymm7, %ymm5
	          vmulpd         none, %ymm5, %ymm5
	          vmovapd        %ymm12, %ymm14
	          vsubpd         UR0, %ymm14, %ymm12
	          vmulpd         none, %ymm12, %ymm12
	          vmovupd        %ymm13, UR0
	          vaddpd         UR1, %ymm7, %ymm13
	          vmovupd        %ymm6, UR1
	          vaddpd         UR1, %ymm15, %ymm6
	          vmovupd        %ymm0, UR3
              vmovupd        UR2, %ymm0
	          vsubpd         UR0, %ymm0, %ymm7
	          vmovupd        UR3, %ymm0
	          vsubpd         UR1, %ymm15, %ymm14
	          vmulpd         none, %ymm14, %ymm14
	          vmovupd        %ymm0, UR3
              vmovupd        UR2, %ymm0
	          vaddpd         UR0, %ymm0, %ymm15
	          vmovupd        UR3, %ymm0
              vmovapd        %ymm2, (%rax)
              vmovapd        %ymm13, (%rax, %rdi, 1)
              vmovapd        %ymm10, (%rbx)
              vmovapd        %ymm5, (%rbx, %rdi, 1)
              vmovapd        %ymm6, (%rcx)
              vmovapd        %ymm9, (%rcx, %rdi, 1)
              vmovapd        %ymm14, (%rdx)
              vmovapd        %ymm1, (%rdx, %rdi, 1)
              pop            %rax
              pop            %rbx
              pop            %rcx
              pop            %rdx
              vmovapd        %ymm0, (%rax)
              vmovapd        %ymm15, (%rax, %rdi, 1)
              vmovapd        %ymm8, (%rbx)
              vmovapd        %ymm7, (%rbx, %rdi, 1)
              vmovapd        %ymm4, (%rcx)
              vmovapd        %ymm11, (%rcx, %rdi, 1)
              vmovapd        %ymm12, (%rdx)
              vmovapd        %ymm3, (%rdx, %rdi, 1)
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
              imul           $7, NLO, %rbx
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
              imul           $3, NLO, %rbx
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
              imul           $5, NLO, %rbx
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
              imul           $1, NLO, %rbx
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
              imul           $6, NLO, %rbx
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
              imul           $2, NLO, %rbx
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
              imul           $4, NLO, %rbx
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

