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

k2rest:       hlt
              mov            X, %rdi
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
skipneg:
              xor            ilo, ilo
k2loop:

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

