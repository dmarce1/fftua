#define       N2o2           %r8
#define       ik             %r9
#define       max_level      %r10
#define       k2             %r11
#define       level          %r12
#define       X              %r13
#define       C              %r14
#define       N              %r15
#define       STACK_SIZE     $0
#define       shft           -8(%rbp)
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
#define       NONE           %ymm14
#define       TWO            %ymm15

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
              bsr            N, %rax
              sub            $4, %rax
              mov            %rax, shft
              mov            N, max_level
              shr            $4, max_level
              bsf            max_level, max_level
              xor            level, level
              xor            k2, k2
L10000:       cmp            k2, N2o2
              jge            L20000
              mov            level, %rcx
              sub            $2, %rcx
              mov            $1, N2o2
              shl            %rcx, N2o2
              shr            N2o2
              mov            k2, ik
              mov            k2, %rbx
              shr            $2, ik
              and            $3, %rbx
              mov            shft, %rcx
              shl            %rcx, %rbx
              or             %rbx, ik
              shl            $2, ik
              nop
L20000:       mov            k2, %rax
              mov            level, %rbx
              inc            %rax
              bsf            %rax, %rcx
              and            $0xfffffffffffffffe, %rcx
              sub            %rcx, %rbx
              shr            %rcx, %rax
              add            $2, level
              shl            $2, k2
              cmp            level, max_level
              cmovl          %rax, k2
              cmovl          %rbx, level
              bsf            k2, %rax
              cmp            level, %rax
              jle            L10000
              pop            %r15
              pop            %r14
              pop            %r13
              pop            %r12
              pop            %rbx
              mov            %rbp, %rsp
              pop            %rbp
              ret
