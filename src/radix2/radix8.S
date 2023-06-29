#include      "common.h"
#define       ER0            %ymm0
#define       ER1            %ymm1
#define       ER2            %ymm2
#define       ER3            %ymm3
#define       EI0            %ymm4
#define       EI1            %ymm5
#define       EI2            %ymm6
#define       EI3            %ymm7
#define       OR0            %ymm8
#define       OR1            %ymm9
#define       OR2            %ymm10
#define       OR3            %ymm11
#define       OI0            %ymm12
#define       OI1            %ymm13
#define       OI2            %ymm14
#define       OI3            %ymm15


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
              shr            $6, %rdi
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
              lea            (%rax, %rsi), %rax
              lea            (%rax, %rsi, 4), %rbx
              lea            (%rax, %rdi, 8), %rcx
              lea            (%rbx, %rdi, 8), %rdx
              vmovapd        (%rax), OR0
              vmovapd        (%rcx), OI0
              vmulpd         SIN1, OI0, ER0
              vmulpd         COS1, OI0, EI0
              vfmsub231pd    COS1, OR0, ER0
              vfmadd231pd    SIN1, OR0, EI0
              vmovapd        ER0, OR0
              vmovapd        EI0, OI0
              vmovapd        (%rax, %rsi, 2), OR1
              vmovapd        (%rcx, %rsi, 2), OI1
              vmulpd         SIN3, OI1, ER0
              vmulpd         COS3, OI1, EI0
              vfmsub231pd    COS3, OR1, ER0
              vfmadd231pd    SIN3, OR1, EI0
              vmovapd        ER0, OR1
              vmovapd        EI0, OI1
              vmovapd        (%rbx), OR2
              vmovapd        (%rdx), OI2
              vmulpd         SIN5, OI2, ER0
              vmulpd         COS5, OI2, EI0
              vfmsub231pd    COS5, OR2, ER0
              vfmadd231pd    SIN5, OR2, EI0
              vmovapd        ER0, OR2
              vmovapd        EI0, OI2
              vmovapd        (%rbx, %rsi, 2), OR3
              vmovapd        (%rdx, %rsi, 2), OI3
              vmulpd         SIN7, OI3, ER0
              vmulpd         COS7, OI3, EI0
              vfmsub231pd    COS7, OR3, ER0
              vfmadd231pd    SIN7, OR3, EI0
              vmovapd        ER0, OR3
              vmovapd        EI0, OI3
              neg            %rsi
              lea            (%rax, %rsi), %rax
              lea            (%rbx, %rsi), %rbx
              lea            (%rcx, %rsi), %rcx
              lea            (%rdx, %rsi), %rdx
              neg            %rsi
              vmovapd        (%rax, %rsi, 2), ER1
              vmovapd        (%rcx, %rsi, 2), EI1
              vmulpd         SIN2, EI1, ER0
              vmulpd         COS2, EI1, EI0
              vfmsub231pd    COS2, ER1, ER0
              vfmadd231pd    SIN2, ER1, EI0
              vmovapd        ER0, ER1
              vmovapd        EI0, EI1
              vmovapd        (%rbx), ER2
              vmovapd        (%rdx), EI2
              vmulpd         SIN4, EI2, ER0
              vmulpd         COS4, EI2, EI0
              vfmsub231pd    COS4, ER2, ER0
              vfmadd231pd    SIN4, ER2, EI0
              vmovapd        ER0, ER2
              vmovapd        EI0, EI2
              vmovapd        (%rbx, %rsi, 2), ER3
              vmovapd        (%rdx, %rsi, 2), EI3
              vmulpd         SIN6, EI3, ER0
              vmulpd         COS6, EI3, EI0
              vfmsub231pd    COS6, ER3, ER0
              vfmadd231pd    SIN6, ER3, EI0
              vmovapd        ER0, ER3
              vmovapd        EI0, EI3
              vmovapd        (%rax), ER0
              vmovapd        (%rcx), EI0
              call           butterfly8
              vmulpd         none, OI0, OI0
              vmulpd         none, OI1, OI1
              vmulpd         none, OI2, OI2
              vmulpd         none, OI3, OI3
              vmovapd        ER0, (%rax)
              vmovapd        ER2, (%rax, %rsi, 2)
              vmovapd        ER1, (%rbx)
              vmovapd        ER3, (%rbx, %rsi, 2)
              vmovapd        OR3, (%rcx)
              vmovapd        OR1, (%rcx, %rsi, 2)
              vmovapd        OR2, (%rdx)
              vmovapd        OR0, (%rdx, %rsi, 2)
              lea            (%rax, %rsi), %rax
              lea            (%rbx, %rsi), %rbx
              lea            (%rcx, %rsi), %rcx
              lea            (%rdx, %rsi), %rdx
              vmovapd        OI0, (%rax)
              vmovapd        OI2, (%rax, %rsi, 2)
              vmovapd        OI1, (%rbx)
              vmovapd        OI3, (%rbx, %rsi, 2)
              vmovapd        EI3, (%rcx)
              vmovapd        EI1, (%rcx, %rsi, 2)
              vmovapd        EI2, (%rdx)
              vmovapd        EI0, (%rdx, %rsi, 2)
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

