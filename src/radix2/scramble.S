#define       X              %r12
#define       N              %r13
#define       N1             %r14
#define       N2             %r15

              .global        fft_scramble

              .text

fft_scramble: push           %r15
              push           %r14
              push           %r13
              push           %r12
              push           %rbx
              mov            %rdi, X
              mov            %rsi, N
              cmp            $4, N
              je             four
              bsf            N, %rcx
              shr            %rcx
              mov            $1, N1
              shl            %rcx, N1
              bsf            N, %rax
              and            $1, %rax
              mov            $1, N2
              cmp            $0, %rax
              mov            $2, %rdx
              cmovne         %rdx, N2
              mov            X, %rdi
              mov            N1, %rsi
              mov            N1, %rdx
              imul           N2, %rdx
              call           scramble_hi
              mov            X, %rdi
              mov            N1, %rsi
              mov            N2, %rdx
              call           transpose_re
              mov            X, %rdi
              mov            N1, %rsi
              mov            N1, %rdx
              imul           N2, %rdx
              call           scramble_hi
              jmp            done
four:         mov            8(X), %rax
              mov            16(X), %rbx
              mov            %rbx, 8(X)
              mov            %rax, 16(X)
done:         pop            %rbx
              pop            %r12
              pop            %r13
              pop            %r14
              pop            %r15
              ret
