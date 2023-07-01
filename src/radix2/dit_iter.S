#define       X              %r8
#define       C              %r9
#define       N              %r10
#define       k2rbits        %r11
#define       k2fbits        %r12
#define       log2NoN1oSS    %r13
#define       NLO            %r14
#define       ilo            %r15
#define       N1             $4
#define       log2N1         $2
#define       SIMD_SIZE      $4
#define       log2SIMD_SIZE  $2
#define       STACK_SIZE     $32
#define       NoN1oSSm1      -8(%rbp)
#define       M              -16(%rbp)
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
#define       COS1           %ymm10
#define       SIN1           %ymm11
#define       COS2           %ymm12
#define       SIN2           %ymm13
#define       tr3            %ymm14
#define       ti3            %ymm15

              .global        dit_iter

              .text

dit_iter:     push           %rbp
              mov            %rsp, %rbp
              sub            STACK_SIZE, %rsp
              push           %rbx
              push           %r12
              push           %r13
              push           %r14
              push           %r15
              mov            %rdi, X
              mov            %rsi, C
              mov            %rdx, N
              bsf            N, log2NoN1oSS
              sub            log2N1, log2NoN1oSS
              sub            log2SIMD_SIZE, log2NoN1oSS
              mov            N, NLO
              shr            log2N1, NLO
              mov            N, %rax
              shr            log2N1, %rax
              shr            log2SIMD_SIZE, %rax
              shl            %rax
              dec            %rax
              mov            %rax, NoN1oSSm1
              mov            $1, k2rbits
              mov            $1, k2fbits
L10000:       mov            k2fbits, %rdi
              bsr            %rdi, %rcx
              mov            $1, %rax
              shl            %rcx, %rax
              not            %rax
              and            %rax, %rdi
              push           %rdi
              mov            %rdi, %rax
              bsr            %rax, %rcx
              mov            $3, %rax
              shl            %rcx, %rax
              dec            %rax
              sub            %rdi, %rax
              mov            %rax, %rsi
              pop            %rdi
              bsr            k2rbits, %rax
              mov            k2rbits, %rdx
              mov            $1, %rax
              shl            %rcx, %rax
              not            %rax
              and            %rax, %rdx
              lea            (C, %rax, 8), %rax
              vbroadcastsd   (%rax), COS1
              vbroadcastsd   (%rax, N, 2), SIN1
              vmulpd         SIN1, SIN1, COS2
              vmulpd         COS1, SIN1, SIN2
              vfmsub231pd    COS1, COS1, COS2
              vfmadd231pd    SIN1, COS1, SIN2
              mov            log2NoN1oSS, %rcx
              add            log2SIMD_SIZE, %rcx
              bsr            k2fbits, %rax
              sub            %rax, %rcx
              mov            $1, %rax
              shl            %rcx, %rax
              mov            %rax, M
              xor            %rcx, %rcx
L15000:       push           %rcx
lea            (X, ilo, 8), %rax
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
              vmovapd        er3, (%rcx)
              vmovapd        er2, (%rax, %rdi, 8)
              vmovapd        ei2, (%rdx)
              vmovapd        ei3, (%rbx, %rdi, 8)
              vfmsub132pd    tr3, er3, er2
              vfnmadd132pd   tr3, ei3, ei2
              pop            %rcx
              add            $4, %rcx
              cmp            %rcx, M
              jne            L15000
              mov            k2fbits, %rdi
              mov            %rdi, %rax
              inc            %rax
              bsf            %rax, %rcx
              and            $0xfffffffffffffffe, %rcx
              shr            %rcx, %rax
              mov            %rdi, %rbx
              shl            log2N1, %rbx
              bsr            %rdi, %rcx
              cmp            %rcx, log2NoN1oSS
              cmovne         %rbx, %rax
              mov            %rax, k2fbits
              mov            k2rbits, %rdi
              mov            %rdi, %rax
              bsr            %rax, %rcx
              mov            $1, %rdx
              shl            %rcx, %rdx
              not            %rdx
              and            %rdx, %rax
              not            %rdx
              dec            %rdx
              mov            %rax, %rbx
              not            %rbx
              and            %rdx, %rbx
              bsr            %rbx, %rcx
              mov            $1, %rbx
              shl            %rcx, %rbx
              dec            %rbx
              mov            %rax, %rcx
              and            %rbx, %rcx
              inc            %rbx
              or             %rbx, %rcx
              mov            %rcx, %rax
              inc            %rdx
              or             %rdx, %rax
              mov            %rdi, %rbx
              shl            log2N1, %rbx
              bsr            %rdi, %rcx
              cmp            %rcx, log2NoN1oSS
              cmovne         %rbx, %rax
              mov            %rax, k2rbits
              cmp            k2fbits, NoN1oSSm1
              jne            L10000
              pop            %r15
              pop            %r14
              pop            %r13
              pop            %r12
              pop            %rbx
              mov            %rbp, %rsp
              pop            %rbp
              ret
