#define       X              %r8
#define       M              %r9
#define       Ix             %r10
#define       Iy             %r11
#define       IMID           %r12
#define       I1             %r13
#define       I2             %r14
#define       N1N2           %r15
#define       M0             $2
#define       STACK_SIZE     $16
#define       N1             -8(%rbp)
#define       N2             -16(%rbp)

              .global        transpose_re_cmplx


              .text

transpose_re_cmplx:
              push           %rbp
              mov            %rsp, %rbp
              sub            STACK_SIZE, %rsp
              push           %r12
              push           %r13
              push           %r14
              push           %r15
              push           %rbx
              mov            %rdi, X
              mov            %rsi, N1
              mov            %rdx, N2
              imul           N1, %rdx
              mov            %rdx, N1N2
              shl            N1N2
              xor            Ix, Ix
              xor            Iy, Iy
              xor            IMID, IMID
              mov            N1, M
              call           next_level
              cmp            $2, N2
              jne            N2one
              xor            Ix, Ix
              xor            Iy, Iy
              mov            $1, IMID
              mov            N1, M
              call           next_level
N2one:        pop            %rbx
              pop            %r15
              pop            %r14
              pop            %r13
              pop            %r12
              mov            %rbp, %rsp
              pop            %rbp
              ret
next_level:   cmp            Ix, Iy
              jl             done
              cmp            M0, M
              jg             next_levels
              cmp            Ix, Iy
              je             on_diagonal
              mov            Ix, %rsi
              mov            Iy, %rdi
              imul           N2, %rsi
              imul           N2, %rdi
              add            IMID, %rsi
              add            IMID, %rdi
              imul           N1, %rsi
              imul           N1, %rdi
              add            Iy, %rsi
              add            Ix, %rdi
              shl            %rsi
              shl            %rdi
              lea            (X, %rsi, 8), %rsi
              lea            (X, %rdi, 8), %rdi
              vmovapd        (%rsi), %ymm0
              vmovapd        (%rsi, N1N2, 8), %ymm1
              vmovapd        (%rdi), %ymm2
              vmovapd        (%rdi, N1N2, 8), %ymm3
              vperm2f128     $0x20, %ymm1, %ymm0, %ymm4
              vperm2f128     $0x31, %ymm1, %ymm0, %ymm1
              vperm2f128     $0x20, %ymm3, %ymm2, %ymm5
              vperm2f128     $0x31, %ymm3, %ymm2, %ymm3
              vmovapd        %ymm4, (%rdi)
              vmovapd        %ymm1, (%rdi, N1N2, 8)
              vmovapd        %ymm5, (%rsi)
              vmovapd        %ymm3, (%rsi, N1N2, 8)
              jmp            done
on_diagonal:  mov            Ix, %rsi
              imul           N2, %rsi
              add            IMID, %rsi
              imul           N1, %rsi
              add            Ix, %rsi
              shl            %rsi
              lea            (X, %rsi, 8), %rsi
              vmovapd        (%rsi), %ymm0
              vmovapd        (%rsi, N1N2, 8), %ymm1
              vperm2f128     $0x20, %ymm1, %ymm0, %ymm2
              vperm2f128     $0x31, %ymm1, %ymm0, %ymm1
              vmovapd        %ymm2, (%rsi)
              vmovapd        %ymm1, (%rsi, N1N2, 8)
              jmp            done
next_levels:  shr            M
              call           next_level
              add            M, Iy
              call           next_level
              sub            M, Iy
              add            M, Ix
              call           next_level
              add            M, Iy
              call           next_level
              sub            M, Iy
              sub            M, Ix
              shl            M
done:         ret





