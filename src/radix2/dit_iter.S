#define       X              %r8
#define       C              %r9
#define       N              %r10
#define       k2rbits        %r11
#define       k2fbits        %r12
#define       log2NoN1       %r13
#define       NLO            %r14
#define       ilo            %r15
#define       N1             $4
#define       log2N1         $2
#define       STACK_SIZE     $32
#define       NoN1m1         -8(%rbp)

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
              bsf            N, log2NoN1
              sub            log2N1, log2NoN1
              mov            N, NLO
              shr            log2N1, NLO
              mov            N, %rax
              shr            log2N1, %rax
              shl            %rax
              dec            %rax
              mov            %rax, NoN1m1
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
              mov            log2NoN1, %rcx
              bsr            k2fbits, %rax
              sub            %rax, %rcx
              mov            $1, %rax
              shl            %rcx, %rax
              mov            %rax, %rcx
L15000:

              loop           L15000
              mov            k2fbits, %rdi
              mov            %rdi, %rax
              inc            %rax
              bsf            %rax, %rcx
              and            $0xfffffffffffffffe, %rcx
              shr            %rcx, %rax
              mov            %rdi, %rbx
              shl            log2N1, %rbx
              bsr            %rdi, %rcx
              cmp            %rcx, log2NoN1
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
              cmp            %rcx, log2NoN1
              cmovne         %rbx, %rax
              mov            %rax, k2rbits
              cmp            k2fbits, NoN1m1
              jne            L10000
              pop            %r15
              pop            %r14
              pop            %r13
              pop            %r12
              pop            %rbx
              mov            %rbp, %rsp
              pop            %rbp
              ret
