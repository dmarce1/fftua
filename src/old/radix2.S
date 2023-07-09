#define       X              %r15
#define       N              %r14
#define       n              %r13
#define       Y              %r12
#define       Nb             $4
#define       STACK_SIZE     $32
#define       N1             -8(%rbp)
#define       N2             -16(%rbp)
#define       stack_pointer  -24(%rbp)


              .global        fft_real_2


              .text

fft_real_2:   push           %rbp
              mov            %rsp, %rbp
              sub            STACK_SIZE, %rsp
              push           %r12
              push           %r13
              push           %r14
              push           %r15
              push           %rbx
              mov            %rdi, X
              mov            %rsi, N
              bsr            N, %rax
              mov            %rax, %rcx
              shr            %rcx
              mov            $1, %rax
              shl            %rcx, %rax
              mov            %rax, N1
              mov            %rax, N2
              and            $1, %rax
              add            $1, %rax
              mov            N2, %rbx
              imul           %rax, %rbx
              mov            %rbx, N2
              mov            N2, %rdx
              imul           Nb, %rdx
              shl            $3, %rdx
              mov            %rsp, stack_pointer
              sub            %rdx, %rsp
              and            $0xffffffffffffffe0, %rsp
              mov            %rsp, Y
              xor            n, n
first_pass:   mov            X, %rdi
              mov            Y, %rsi
              mov            N1, %rdx
              mov            N2, %rcx
              mov            n, %r8
              call           transpose_in
              lea            (Y), %rdi
              mov            N2, %rsi
              xor            %rdx, %rdx
              call           dit_nr_recur
              mov            N2, %rax
              lea            (Y, %rax, 8), %rdi
              mov            N2, %rsi
              xor            %rdx, %rdx
              call           dit_nr_recur
              imul           $2, N2, %rax
              lea            (Y, %rax, 8), %rdi
              mov            N2, %rsi
              xor            %rdx, %rdx
              call           dit_nr_recur
              imul           $3, N2, %rax
              lea            (Y, %rax, 8), %rdi
              mov            N2, %rsi
              xor            %rdx, %rdx
              call           dit_nr_recur
              mov            Y, %rdi
              mov            X, %rsi
              mov            N1, %rdx
              mov            N2, %rcx
              mov            n, %r8
              call           transpose_out
              add            Nb, n
              cmp            n, N1
              jne            first_pass
              xor            n, n
second_pass:  mov            X, %rdi
              mov            N, %rsi
              bsf            N2, %rdx
              call           dit_nr_recur
              mov            X, %rdi
              mov            N, %rsi
              call           scramble
              mov            stack_pointer, %rsp
              pop            %rbx
              pop            %r15
              pop            %r14
              pop            %r13
              pop            %r12
              mov            %rbp, %rsp
              pop            %rbp
              ret
