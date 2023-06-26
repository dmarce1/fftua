#include      "registers.h"

              .global        addr_adjust

              .text

addr_adjust:  mov            N2, %rdx
              dec            %rdx
              mov            X, %rax
              bsf            NLO, %rcx
              sub            X0, %rax
              add            $3, %rcx
              add            %rdi, %rcx
              shr            %rcx, %rax
              mov            $1, %rbx
              mov            %rax, %rsi
              bsf            N2, %rcx
              not            %rax
              dec            %rcx
              and            %rdx, %rax
              shl            %rcx, %rbx
rev_add_b:    test           %rbx, %rax
              jz             rev_add_e
              xor            %rbx, %rax
              shr            %rbx
              jmp            rev_add_b
rev_add_e:    or             %rbx, %rax
              sub            %rsi, %rax
              jl             skip
              imul           NLO, %rax
              mov            %rdi, %rcx
              shl            %rcx, %rax
              mov            N2, %rsi
              shr            %rsi
              cmp            k2, %rsi
              jg             done
              lea            (X, %rax, 8), X
              neg            %rax
              mov            N2, %rsi
              sub            k2, %rsi
              mov            %rsi, k2
              jmp            done
skip:         xor            %rax, %rax
              not            %rax
done:         ret
