#define       X              %r8
#define       C              %r9
#define       N              %r10
#define       N0             %r11
#define       NLO            %r12
#define       k2             %r13
#define       ilo            %r14
#define       N2             %r15
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
#define       STACK_SIZE     $32
#define       X0             -8(%rbp)


              .global        dit_nr_entry


              .text

dit_nr_entry: push           %rbp
              mov            %rsp, %rbp
              sub            STACK_SIZE, %rsp
              push           %rbx
              push           %r12
              push           %r13
              push           %r14
              push           %r15
              mov            %rdi, X0
              mov            %rdx, N0
              mov            %rdx, N
              mov            %rsi, C
              mov            %rdi, X
              mov            $1, N2
              xor            k2, k2
              call           dit_nr
              pop            %r15
              pop            %r14
              pop            %r13
              pop            %r12
              pop            %rbx
              mov            %rbp, %rsp
              pop            %rbp
              ret
dit_nr:       bsf            N, %rax
              test           $1, %rax
              jnz            radix8
              mov            N, NLO
              shr            $2, NLO
              cmp            $0, k2
              je             keq0
              jmp            k2rest
keq0:         xor            ilo, ilo
              mov            NLO, %rsi
              shl            $3, %rsi
              cmp            $1, N2
              jle            n21
              vmovapd        tw45, ti3
              vmovapd        none, ti2
k0loop:       lea            (X, ilo, 8), %rax
              lea            (%rax, %rsi, 2), %rbx
              lea            (X, %rsi, 4), %rcx
              lea            (%rcx, ilo, 8), %rcx
              lea            (%rcx, %rsi, 2), %rdx
              vmovapd        (%rax), er0
              vmovapd        (%rbx), er2
              vmovapd        (%rax, %rsi), er1
              vmovapd        (%rbx, %rsi), er3
              vmovapd        (%rcx), ei0
              vmovapd        (%rdx), ei2
              vmovapd        (%rcx, %rsi), ei1
              vmovapd        (%rdx, %rsi), ei3
              vaddpd         ei3, ei1, tr0
              vsubpd         ei3, ei1, tr2
              vmulpd         ti3, tr2, tr1
              vmulpd         ti3, tr0, tr3
              vmovapd        ei0, tr0
              vmovapd        ei2, tr2
              vaddpd         tr1, tr0, ei0
              vaddpd         tr3, tr2, ei3
              vsubpd         tr1, tr0, ei1
              vsubpd         tr3, tr2, ei2
              vmulpd         ti2, ei3, ei3
              vaddpd         er2, er0, tr0
              vsubpd         er2, er0, tr2
              vaddpd         er3, er1, tr1
              vsubpd         er1, er3, tr3
              vaddpd         tr1, tr0, er0
              vsubpd         tr1, tr0, er2
              vmovapd        ei0, (%rcx)
              vmovapd        ei1, (%rdx)
              vmovapd        ei2, (%rcx, %rsi)
              vmovapd        ei3, (%rdx, %rsi)
              vmovapd        er0, (%rax)
              vmovapd        er2, (%rax, %rsi)
              vmovapd        tr2, (%rbx)
              vmovapd        tr3, (%rbx, %rsi)
              add            $4, ilo
              cmp            ilo, NLO
              jg             k0loop
              jmp            rcalls
n21:          mov            $4, %rdi
n21loop:      lea            (X, ilo, 8), %rax
              lea            (%rax, %rdi, 8), %rcx
              lea            (%rax, %rsi, 2), %rbx
              lea            (%rcx, %rsi, 2), %rdx
              vmovapd        (%rax), er0
              vmovapd        (%rbx), er2
              vmovapd        (%rax, %rsi), er1
              vmovapd        (%rbx, %rsi), er3
              vmovapd        (%rcx), ei0
              vmovapd        (%rdx), ei2
              vmovapd        (%rcx, %rsi), ei1
              vmovapd        (%rdx, %rsi), ei3
              vaddpd         er2, er0, tr0
              vsubpd         er2, er0, tr2
              vaddpd         er3, er1, tr1
              vsubpd         er1, er3, tr3
              vaddpd         tr1, tr0, er0
              vsubpd         tr1, tr0, er2
              vaddpd         ei2, ei0, ti0
              vsubpd         ei2, ei0, ti2
              vaddpd         ei3, ei1, ti1
              vsubpd         ei1, ei3, ti3
              vaddpd         ti1, ti0, ei0
              vsubpd         ti1, ti0, ei2
              vmovapd        er0, (%rax)
              vmovapd        er2, (%rax, %rsi)
              vmovapd        tr2, (%rbx)
              vmovapd        tr3, (%rbx, %rsi)
              vmovapd        ei0, (%rcx)
              vmovapd        ei2, (%rcx, %rsi)
              vmovapd        ti2, (%rdx)
              vmovapd        ti3, (%rdx, %rsi)
              add            $8, ilo
              cmp            ilo, NLO
              jg             n21loop
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
              jle            rcalls
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
skipneg:      mov            k2, %rax
              imul           NLO, %rax
              lea            (C, %rax, 8), %rax
              vbroadcastsd   (%rax), tr1
              vbroadcastsd   (%rax, N0, 2), ti1
              vmulpd         ti1, ti1, tr2
              vmulpd         tr1, ti1, ti2
              vfmsub231pd    tr1, tr1, tr2
              vfmadd231pd    ti1, tr1, ti2
              vmovapd        two, tr3
              xor            ilo, ilo
k2loop:       lea            (X, ilo, 8), %rax
              vmovapd        (%rax), er0
              vmovapd        (%rax, %rdi, 8), ei0
              lea            (%rax, NLO, 8), %rbx
              lea            (%rbx, NLO, 8), %rcx
              vmovapd        (%rcx), er2
              vmovapd        er0, tr0
              vmovapd        (%rcx, %rdi, 8), ei2
              vmovapd        ei0, ti0
              lea            (%rcx, NLO, 8), %rdx
              vfmadd231pd    ti2, ei2, tr0
              vmovapd        (%rbx), er1
              vfnmadd231pd   ti2, er2, ti0
              vmovapd        (%rdx), er3
              vfnmadd132pd   tr2, tr0, er2
              vmovapd        (%rbx, %rdi, 8), ei1
              vfnmadd132pd   tr2, ti0, ei2
              vmovapd        (%rdx, %rdi, 8), ei3
              vfmsub132pd    tr3, er2, er0
              vfmsub132pd    tr3, ei2, ei0
              vmovapd        er1, tr0
              vmovapd        ei1, ti0
              vfmadd231pd    ti2, ei3, tr0
              vfnmadd231pd   ti2, er3, ti0
              vfnmadd132pd   tr2, tr0, er3
              vfnmadd132pd   tr2, ti0, ei3
              vfmsub132pd    tr3, er3, er1
              vfmsub132pd    tr3, ei3, ei1
              vmovapd        er0, tr0
              vmovapd        ei0, ti0
              vfmadd231pd    ti1, ei1, tr0
              vfnmadd231pd   ti1, er1, ti0
              vfnmadd132pd   tr1, tr0, er1
              vfmsub132pd    tr1, ti0, ei1
              vfmsub132pd    tr3, er1, er0
              vfmadd132pd    tr3, ei1, ei0
              vmovapd        er2, tr0
              vmovapd        ei2, ti0
              vfmadd231pd    tr1, ei3, tr0
              vmovapd        er0, (%rax)
              vfnmadd231pd   tr1, er3, ti0
              vmovapd        ei0, (%rdx, %rdi, 8)
              vfmadd132pd    ti1, tr0, er3
              vmovapd        er1, (%rcx, %rdi, 8)
              vfmadd132pd    ti1, ti0, ei3
              vmovapd        ei1, (%rbx)
              vfmsub132pd    tr3, er3, er2
              vmovapd        er3, (%rcx)
              vfnmadd132pd   tr3, ei3, ei2
              vmovapd        er2, (%rax, %rdi, 8)
              vmovapd        ei2, (%rdx)
              vmovapd        ei3, (%rbx, %rdi, 8)
              add            $4, ilo
              cmp            ilo, NLO
              jg             k2loop
              pop            k2
              pop            X
rcalls:       cmp            $4, NLO
              jle            return
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
              call           dit_nr
              pop            k2
              pop            X
              call           dit_nr
              pop            k2
              pop            X
              call           dit_nr
              pop            k2
              pop            X
              call           dit_nr
              pop            k2
              pop            X
              mov            N, %rax
              shl            $2, %rax
              mov            %rax, N
              shr            $2, N2
              jmp            return
radix8:       mov            N, NLO
              shr            $3, NLO
              xor            ilo, ilo
              vmovapd        none, %ymm15
              vmovapd        tw45, %ymm14
              mov            NLO, %rdi
              shl            $3, %rdi
radix8lp:     lea            (X, ilo, 8), %rax
              lea            (%rax, %rdi, 2), %rbx
              lea            (%rbx, %rdi, 2), %rcx
              lea            (%rcx, %rdi, 2), %rdx
              vmovapd        (%rax), %ymm0
              vmovapd        (%rax, %rdi), %ymm1
              vmovapd        (%rbx), %ymm2
              vmovapd        (%rbx, %rdi), %ymm3
              vmovapd        (%rcx), %ymm4
              vmovapd        (%rcx, %rdi), %ymm5
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
	          vmulpd         %ymm15, %ymm3, %ymm6
	          vmovapd        %ymm4, %ymm3
	          vsubpd         %ymm8, %ymm5, %ymm4
	          vmulpd         %ymm14, %ymm7, %ymm7
	          vmulpd         %ymm14, %ymm1, %ymm1
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
              vmovapd        %ymm4, (%rax, %rdi)
              vmovapd        %ymm2, (%rbx)
              vmovapd        %ymm6, (%rbx, %rdi)
              vmovapd        %ymm1, (%rcx)
              vmovapd        %ymm5, (%rcx, %rdi)
              vmovapd        %ymm3, (%rdx)
              vmovapd        %ymm7, (%rdx, %rdi)
              add            $4, ilo
              cmp            ilo, NLO
              jg             radix8lp
              cmp            $4, NLO
              jle            return
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
              call           dit_nr
              pop            k2
              pop            X
              call           dit_nr
              pop            k2
              pop            X
              call           dit_nr
              pop            k2
              pop            X
              call           dit_nr
              pop            k2
              pop            X
              call           dit_nr
              pop            k2
              pop            X
              call           dit_nr
              pop            k2
              pop            X
              call           dit_nr
              pop            k2
              pop            X
              call           dit_nr
              pop            k2
              pop            X
              mov            N, %rax
              shl            $3, %rax
              mov            %rax, N
              shr            $3, N2
return:       ret

