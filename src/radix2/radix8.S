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
#define       TMP0           (%rsp)
#define       TMP1           32(%rsp)
#define       TMP2           64(%rsp)
#define       TMP3           96(%rsp)
#define       TMP4           128(%rsp)
#define       IDIST          160(%rsp)
#define       COS1           192(%rsp)
#define       COS2           200(%rsp)
#define       COS3           208(%rsp)
#define       COS4           216(%rsp)
#define       COS5           224(%rsp)
#define       COS6           232(%rsp)
#define       COS7           240(%rsp)
#define       COS8           248(%rsp)
#define       SIN1           256(%rsp)
#define       SIN2           264(%rsp)
#define       SIN3           272(%rsp)
#define       SIN4           280(%rsp)
#define       SIN5           288(%rsp)
#define       SIN6           296(%rsp)
#define       SIN7           304(%rsp)
#define       SIN8           312(%rsp)
#define       STACK_SIZE     $320

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
signbit:      .quad          0x8000000000000000
              .quad          0x8000000000000000
              .quad          0x8000000000000000
              .quad          0x8000000000000000


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
	          vxorpd         signbit,%ymm3, %ymm6
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
	          vxorpd         signbit,%ymm6, %ymm7
	          vmovapd        %ymm3, %ymm6
	          vsubpd         %ymm8, %ymm5, %ymm3
	          vsubpd         %ymm10, %ymm4, %ymm4
	          vmovapd        %ymm2, %ymm8
	          vaddpd         %ymm1, %ymm6, %ymm2
	          vsubpd         %ymm9, %ymm8, %ymm5
	          vmovapd        %ymm1, %ymm10
	          vsubpd         %ymm6, %ymm10, %ymm1
	          vaddpd         %ymm8, %ymm9, %ymm9
	          vxorpd         signbit,%ymm9, %ymm6
	          vmovapd        %ymm4, %ymm8
	          vxorpd         signbit,%ymm8, %ymm4
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
skipneg:      push           %rbp
              mov            %rsp, %rbp
              sub            STACK_SIZE, %rsp
              and            $0xffffffffffffffc0, %rsp
              mov            k2, %rax
              imul           NLO, %rax
              lea            (C, %rax, 8), %rax
              vmovq          (%rax), %xmm0
              vmovq          (%rax, N0, 2), %xmm1
              vmulsd         %xmm1, %xmm1, %xmm2
              vmulsd         %xmm0, %xmm1, %xmm3
              vfmsub231sd    %xmm0, %xmm0, %xmm2
              vfmadd231sd    %xmm1, %xmm0, %xmm3
              vbroadcastsd   %xmm2, %ymm2
              vbroadcastsd   %xmm3, %ymm3
              vbroadcastsd   %xmm0, %ymm0
              vbroadcastsd   %xmm1, %ymm1
              vmovq          %xmm2, COS1
              vmovq          %xmm3, SIN1
              vmovhpd        COS1, %xmm0, %xmm0
              vmovhpd        SIN1, %xmm1, %xmm1
              vmulpd         %xmm1, %xmm3, %xmm4
              vmulpd         %xmm0, %xmm3, %xmm5
              vfmsub231pd    %xmm0, %xmm2, %xmm4
              vfmadd231pd    %xmm1, %xmm2, %xmm5
              vinsertf128    $1, %xmm4, %ymm0, %ymm2
              vinsertf128    $1, %xmm5, %ymm1, %ymm3
              vpermpd        $255, %ymm2, %ymm0
              vpermpd        $255, %ymm3, %ymm1
              vmulpd         %ymm1, %ymm3, %ymm4
              vmulpd         %ymm0, %ymm3, %ymm5
              vfmsub231pd    %ymm0, %ymm2, %ymm4
              vfmadd231pd    %ymm1, %ymm2, %ymm5
              vmovapd        %ymm2, COS1
              vmovapd        %ymm3, SIN1
              vmovapd        %ymm4, COS5
              vmovapd        %ymm5, SIN5
              mov            %rdi, IDIST
              xor            ilo, ilo
k2loop:       imul           $8, NLO, %rax
              imul           $24, NLO, %rbx
              imul           $40, NLO, %rcx
              imul           $56, NLO, %rdx
              lea            (X, ilo, 8), %rsi
              mov            IDIST, %rdi
              lea            (%rsi, %rdi, 8), %rdi
              vmovapd        (%rsi, %rax, 1), %ymm0
              vmovapd        (%rdi, %rax, 1), %ymm1
              vbroadcastsd   COS1, %ymm14
              vbroadcastsd   SIN1, %ymm15
              vmulpd         %ymm15, %ymm1, %ymm2
              vmulpd         %ymm14, %ymm1, %ymm3
              vfmsub231pd    %ymm14, %ymm0, %ymm2
              vfmadd231pd    %ymm15, %ymm0, %ymm3
              vmovapd        (%rsi, %rax, 2), %ymm0
              vmovapd        (%rdi, %rax, 2), %ymm1
              vbroadcastsd   COS2, %ymm14
              vbroadcastsd   SIN2, %ymm15
              vmulpd         %ymm15, %ymm1, %ymm4
              vmulpd         %ymm14, %ymm1, %ymm5
              vfmsub231pd    %ymm14, %ymm0, %ymm4
              vfmadd231pd    %ymm15, %ymm0, %ymm5
              vmovapd        (%rsi, %rbx, 1), %ymm0
              vmovapd        (%rdi, %rbx, 1), %ymm1
              vbroadcastsd   COS3, %ymm14
              vbroadcastsd   SIN3, %ymm15
              vmulpd         %ymm15, %ymm1, %ymm6
              vmulpd         %ymm14, %ymm1, %ymm7
              vfmsub231pd    %ymm14, %ymm0, %ymm6
              vfmadd231pd    %ymm15, %ymm0, %ymm7
              vmovapd        (%rsi, %rax, 4), %ymm0
              vmovapd        (%rdi, %rax, 4), %ymm1
              vbroadcastsd   COS4, %ymm14
              vbroadcastsd   SIN4, %ymm15
              vmulpd         %ymm15, %ymm1, %ymm8
              vmulpd         %ymm14, %ymm1, %ymm9
              vfmsub231pd    %ymm14, %ymm0, %ymm8
              vfmadd231pd    %ymm15, %ymm0, %ymm9
              vmovapd        (%rsi, %rcx, 1), %ymm0
              vmovapd        (%rdi, %rcx, 1), %ymm1
              vbroadcastsd   COS5, %ymm14
              vbroadcastsd   SIN5, %ymm15
              vmulpd         %ymm15, %ymm1, %ymm10
              vmulpd         %ymm14, %ymm1, %ymm11
              vfmsub231pd    %ymm14, %ymm0, %ymm10
              vfmadd231pd    %ymm15, %ymm0, %ymm11
              vmovapd        (%rsi, %rbx, 2), %ymm0
              vmovapd        (%rdi, %rbx, 2), %ymm1
              vbroadcastsd   COS6, %ymm14
              vbroadcastsd   SIN6, %ymm15
              vmulpd         %ymm15, %ymm1, %ymm12
              vmulpd         %ymm14, %ymm1, %ymm13
              vfmsub231pd    %ymm14, %ymm0, %ymm12
              vfmadd231pd    %ymm15, %ymm0, %ymm13
              vmovapd        %ymm12, TMP0
              vmovapd        %ymm13, TMP1
              vmovapd        (%rsi, %rdx, 1), %ymm0
              vmovapd        (%rdi, %rdx, 1), %ymm1
              vbroadcastsd   COS7, %ymm12
              vbroadcastsd   SIN7, %ymm13
              vmulpd         %ymm13, %ymm1, %ymm14
              vmulpd         %ymm12, %ymm1, %ymm15
              vfmsub231pd    %ymm12, %ymm0, %ymm14
              vfmadd231pd    %ymm13, %ymm0, %ymm15
              vmovapd        TMP0, %ymm12
              vmovapd        TMP1, %ymm13
              vmovapd        (%rsi), %ymm0
              vmovapd        (%rdi), %ymm1
  	          vmovapd        %ymm8, TMP3
              vmovapd        %ymm9, TMP4
	          vaddpd         %ymm0, %ymm8, %ymm8
	          vaddpd         %ymm1, %ymm9, %ymm9
	          vmovapd        %ymm8, TMP0
	          vmovapd        %ymm9, TMP1
	          vmovapd        TMP3, %ymm8
	          vmovapd        TMP4, %ymm9
 	          vsubpd         %ymm8, %ymm0, %ymm0
              vsubpd         %ymm9, %ymm1, %ymm1
	          vaddpd         %ymm4, %ymm12, %ymm8
	          vaddpd         %ymm5, %ymm13, %ymm9
	          vsubpd         %ymm12, %ymm4, %ymm4
	          vsubpd         %ymm13, %ymm5, %ymm5
	          vaddpd         TMP0, %ymm8, %ymm12
	          vaddpd         TMP1, %ymm9, %ymm13
	          vmovapd        %ymm0, TMP3
              vmovapd        %ymm1, TMP4
              vmovapd        TMP0, %ymm0
              vmovapd        TMP1, %ymm1
              vsubpd         %ymm8, %ymm0, %ymm0
              vsubpd         %ymm9, %ymm1, %ymm1
              vmovapd        %ymm0, TMP0
              vmovapd        %ymm1, TMP1
	          vmovapd        TMP3, %ymm0
	          vmovapd        TMP4, %ymm1
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
	          vmulpd         C2, %ymm2, %ymm2
	          vmulpd         C2, %ymm6, %ymm6
	          vsubpd         %ymm2, %ymm6, %ymm10
	          vaddpd         %ymm6, %ymm2, %ymm2
	          vmovapd        %ymm0, %ymm6
	          vaddpd         %ymm14, %ymm12, %ymm0
	          vmovapd        %ymm1, %ymm7
	          vaddpd         %ymm15, %ymm13, %ymm1
              vmovapd        %ymm1, (%rdi, %rdx, 1)
	          vmovapd        %ymm8, %ymm1
	          vsubpd         %ymm14, %ymm12, %ymm8
              vmovapd        %ymm8, (%rdi, %rbx, 2)
	          vmovapd        %ymm9, %ymm12
	          vsubpd         %ymm15, %ymm13, %ymm9
	          vmovapd        %ymm2, %ymm13
	          vaddpd         %ymm1, %ymm3, %ymm2
              vmovapd        %ymm2, (%rsi, %rax, 4)
	          vmovapd        %ymm3, %ymm14
	          vaddpd         %ymm11, %ymm12, %ymm3
              vmovapd        %ymm3, (%rdi, %rbx, 1)
	          vmovapd        %ymm10, %ymm15
	          vsubpd         %ymm1, %ymm14, %ymm10
	          vmovapd        signbit, %ymm8
	          vxorpd         %ymm8, %ymm10, %ymm10
              vmovapd        %ymm10, (%rdi, %rax, 2)
	          vmovapd        %ymm11, %ymm1
	          vsubpd         %ymm1, %ymm12, %ymm11
	          vmovapd        %ymm4, %ymm1
	          vmovapd        TMP0, %ymm2
	          vaddpd         %ymm2, %ymm5, %ymm4
              vmovapd        %ymm4, (%rsi, %rax, 2)
	          vmovapd        %ymm5, %ymm12
	          vmovapd        TMP1, %ymm3
	          vsubpd         %ymm3, %ymm7, %ymm5
	          vxorpd         %ymm8, %ymm5, %ymm5
              vmovapd        %ymm5, (%rdi, %rcx, 1)
	          vmovapd        %ymm12, %ymm14
	          vsubpd         %ymm2, %ymm14, %ymm12
	          vxorpd         %ymm8, %ymm12, %ymm12
              vmovapd        %ymm12, (%rdi, %rax, 4)
	          vmovapd        %ymm13, %ymm2
	          vaddpd         %ymm3, %ymm7, %ymm13
	          vmovapd        %ymm6, %ymm3
	          vaddpd         %ymm3, %ymm15, %ymm6
              vmovapd        %ymm6, (%rsi, %rbx, 2)
	          vmovapd        %ymm0, %ymm10
              vmovapd        %ymm1, %ymm0
	          vsubpd         %ymm2, %ymm0, %ymm7
              vmovapd        %ymm7, (%rdi, %rax, 1)
	          vsubpd         %ymm3, %ymm15, %ymm14
	          vxorpd         %ymm8, %ymm14, %ymm14
              vmovapd        %ymm14, (%rdi)
              vmovapd        %ymm1, %ymm0
	          vaddpd         %ymm2, %ymm0, %ymm15
	          vmovapd        %ymm10, %ymm0
              vmovapd        %ymm0, (%rsi)
              vxorpd         %ymm8, %ymm9, %ymm9
              vmovapd        %ymm9, (%rsi, %rax, 1)
              vxorpd         %ymm8, %ymm11, %ymm11
              vmovapd        %ymm11, (%rsi, %rcx, 1)
              vxorpd         %ymm8, %ymm13, %ymm13
              vmovapd        %ymm13, (%rsi, %rbx, 1)
              vxorpd         %ymm8, %ymm15, %ymm15
              vmovapd        %ymm15, (%rsi, %rdx, 1)
              add            $4, ilo
              cmp            ilo, NLO
              jg             k2loop
              mov            %rbp, %rsp
              pop            %rbp
              pop            k2
              pop            X
rcalls:       cmp            $4, NLO
              jle            done1
              push           X
              push           k2
              bsf            N2, %rcx
              mov            $7, %rax
              imul           $7, NLO, %rbx
              shl            %rcx, %rax
              or             k2, %rax
              lea            (X, %rbx, 8), %rcx
              push           %rcx
              push           %rax
              bsf            N2, %rcx
              mov            $6, %rax
              imul           $3, NLO, %rbx
              shl            %rcx, %rax
              or             k2, %rax
              lea            (X, %rbx, 8), %rcx
              push           %rcx
              push           %rax
              bsf            N2, %rcx
              mov            $5, %rax
              imul           $5, NLO, %rbx
              shl            %rcx, %rax
              or             k2, %rax
              lea            (X, %rbx, 8), %rcx
              push           %rcx
              push           %rax
              bsf            N2, %rcx
              mov            $4, %rax
              imul           $1, NLO, %rbx
              shl            %rcx, %rax
              or             k2, %rax
              lea            (X, %rbx, 8), %rcx
              push           %rcx
              push           %rax
              bsf            N2, %rcx
              mov            $3, %rax
              imul           $6, NLO, %rbx
              shl            %rcx, %rax
              or             k2, %rax
              lea            (X, %rbx, 8), %rcx
              push           %rcx
              push           %rax
              bsf            N2, %rcx
              mov            $2, %rax
              imul           $2, NLO, %rbx
              shl            %rcx, %rax
              or             k2, %rax
              lea            (X, %rbx, 8), %rcx
              push           %rcx
              push           %rax
              bsf            N2, %rcx
              mov            $1, %rax
              imul           $4, NLO, %rbx
              shl            %rcx, %rax
              or             k2, %rax
              lea            (X, %rbx, 8), %rcx
              push           %rcx
              push           %rax
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

