#include      "registers.h"
#define       RADIX          $2
#define       log2RADIX      $1


              .text

radix2:       mov            N, N2
              shr            log2RADIX, N2
              cmp            $0, k2
              je             k2eq0
              mov            N2, %rax
              shr            %rax
              cmp            %rax, k2
              je             k2Ny
              jmp            k2rest
k2eq0:        xor            ilo, ilo
k2eq0_loop:   lea            (X, ilo, 8), %rax
              lea            (%rax, NLO, 4), %rbx
              vmovapd        (%rax), er0
              vmovapd        (%rbx), er1
              vaddpd         er1, er0, tr0
              vsubpd         er1, er0, tr1
              vmovapd        er0, (%rax)
              vmovapd        er1, (%rbx)
              add            SIMD_SIZE, ilo
              cmp            ilo, NLO
              jne            k2eq0_loop
              jmp            finish
k2Ny:         xor            ilo, ilo
k2Ny_loop:    lea            (X, ilo, 8), %rax
              lea            (%rax, NLO, 4), %rbx
              vmovapd        (%rbx), er1
              vmulpd         none, er1, er1
              vmovapd        er1, (%rbx)
              add            SIMD_SIZE, ilo
              cmp            ilo, NLO
              jne            k2Ny_loop
              jmp            finish
k2rest:       push           X
              push           k2
              mov            log2RADIX, %rdi
              call           addr_adjust
              cmp            $0xffffffffffffffff, %rax
              je             skip_k2rest
              mov            %rax, %rdi
		      mov            k2, %rax
              imul           NLO, %rax
              lea            (C, %rax, 8), %rax
              mov            N0, %rdx
              vbroadcastsd   (%rax), tr1
              vbroadcastsd   (%rax, %rdx, 2), ti1
              xor            ilo, ilo
k2rest_loop:  lea            (X, ilo, 8), %rax
              lea            (%rax, NLO, 4), %rbx
              vmovapd        (%rax), er0
              vmovapd        (%rbx), er1
              vmovapd        (%rax, %rdi, 8), ei0
              vmovapd        (%rbx, %rdi, 8), ei1
              vmulpd         ti1, ei1, tr0
              vmulpd         tr1, ei1, ti0
              vfmsub231pd    tr1, er1, tr0
              vfmadd231pd    ti1, er1, ti0
              vmovapd        tr0, er1
              vmovapd        ti0, ei1
              vaddpd         er1, er0, tr0
              vsubpd         er1, er0, tr1
              vaddpd         ei1, ei0, ti0
              vsubpd         ei1, ei0, ti1
              vmovapd        tr0, (%rax)
              vmovapd        ti0, (%rbx)
              vmovapd        tr1, (%rax, %rdi, 8)
              vmovapd        ti1, (%rbx, %rdi, 8)
              add            SIMD_SIZE, ilo
              cmp            ilo, NLO
              jne            k2rest_loop
skip_k2rest:  pop            k2
              pop            X
finish:       push           X
              push           k2
              bsf            N2, %rcx
              mov            $1, %rax
              shl            %rcx, %rax
              or             k2, %rax
              lea            (X, NLO, 4), %rbx
              push           %rbx
              push           %rax
              shl            log2RADIX, N
              shr            log2RADIX, NLO
              call           radix2
              pop            k2
              pop            X
              call           radix2
              pop            k2
              pop            X
              shl            log2RADIX, NLO
              shr            log2RADIX, N
skip_finish:  ret







