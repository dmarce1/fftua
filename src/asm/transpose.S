#define       X              %r8
#define       M              %r9
#define       Ix             %r10
#define       Iy             %r11
#define       IMID           %r12
#define       I1             %r13
#define       I2             %r14
#define       N1N2           %r15
#define       M0             $4
#define       STACK_SIZE     $16
#define       N1             -8(%rbp)
#define       N2             -16(%rbp)

              .global        transpose_complex


              .text

transpose_complex:
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
              xor            Ix, Ix
              xor            Iy, Iy
              mov            N1, M
              call           next_level
              pop            %rbx
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
              imul           $2, N1N2, %rax
              imul           $3, N1N2, %rbx
              xor            IMID, IMID
              cmp            Ix, Iy
              je             on_diagonal
main_loop:    mov            Ix, %rsi
              mov            Iy, %rdi
              imul           N2, %rsi
              imul           N2, %rdi
              add            IMID, %rsi
              add            IMID, %rdi
              imul           N1, %rsi
              imul           N1, %rdi
              add            Iy, %rsi
              add            Ix, %rdi
              lea            (X, %rsi, 8), %rsi
              lea            (X, %rdi, 8), %rdi
              vmovapd        (%rsi), %ymm0
              vmovapd        (%rsi, N1N2, 8), %ymm1
              vmovapd        (%rdi), %ymm2
              vmovapd        (%rdi, N1N2, 8), %ymm3
              vpermpd        $78, %ymm0, %ymm4
              vinsertf128    $1, %xmm1, %ymm0, %ymm0
              vinsertf128    $0, %xmm4, %ymm1, %ymm1
              vpermpd        $78, %ymm2, %ymm4
              vinsertf128    $1, %xmm3, %ymm2, %ymm2
              vinsertf128    $0, %xmm4, %ymm3, %ymm3
              vmovapd        %ymm0, (%rdi)
              vmovapd        %ymm1, (%rdi, N1N2, 8)
              vmovapd        %ymm2, (%rsi)
              vmovapd        %ymm3, (%rsi, N1N2, 8)
              vmovapd        32(%rsi, %rax, 8), %ymm0
              vmovapd        32(%rsi, %rbx, 8), %ymm1
              vmovapd        32(%rdi, %rax, 8), %ymm2
              vmovapd        32(%rdi, %rbx, 8), %ymm3
              vpermpd        $78, %ymm0, %ymm4
              vinsertf128    $1, %xmm1, %ymm0, %ymm0
              vinsertf128    $0, %xmm4, %ymm1, %ymm1
              vpermpd        $78, %ymm2, %ymm4
              vinsertf128    $1, %xmm3, %ymm2, %ymm2
              vinsertf128    $0, %xmm4, %ymm3, %ymm3
              vmovapd        %ymm0, 32(%rdi, %rax, 8)
              vmovapd        %ymm1, 32(%rdi, %rbx, 8)
              vmovapd        %ymm2, 32(%rsi, %rax, 8)
              vmovapd        %ymm3, 32(%rsi, %rbx, 8)
              vmovapd        32(%rsi), %ymm0
              vmovapd        32(%rsi, N1N2, 8), %ymm1
              vmovapd        (%rsi, %rax, 8), %ymm2
              vmovapd        (%rsi, %rbx, 8), %ymm3
              vmovapd        32(%rdi), %ymm4
              vmovapd        32(%rdi, N1N2, 8), %ymm5
              vmovapd        (%rdi, %rax, 8), %ymm6
              vmovapd        (%rdi, %rbx, 8), %ymm7
              vpermpd        $78, %ymm0, %ymm8
              vinsertf128    $1, %xmm1, %ymm0, %ymm0
              vinsertf128    $0, %xmm8, %ymm1, %ymm1
              vpermpd        $78, %ymm2, %ymm8
              vinsertf128    $1, %xmm3, %ymm2, %ymm2
              vinsertf128    $0, %xmm8, %ymm3, %ymm3
              vpermpd        $78, %ymm0, %ymm8
              vinsertf128    $1, %xmm5, %ymm4, %ymm4
              vinsertf128    $0, %xmm8, %ymm5, %ymm5
              vpermpd        $78, %ymm6, %ymm8
              vinsertf128    $1, %xmm7, %ymm6, %ymm6
              vinsertf128    $0, %xmm8, %ymm7, %ymm7
              vmovapd        %ymm3, 32(%rdi, %rax, 8)
              vmovapd        %ymm2, 32(%rdi, %rbx, 8)
              vmovapd        32(%rdi), %ymm1
              vmovapd        32(%rdi, N1N2, 8), %ymm0
              vmovapd        %ymm7, 32(%rsi, %rax, 8)
              vmovapd        %ymm6, 32(%rsi, %rbx, 8)
              vmovapd        32(%rsi), %ymm5
              vmovapd        32(%rsi, N1N2, 8), %ymm4
              jmp            skip_diagonal
on_diagonal:  mov            Ix, %rsi
              imul           N2, %rsi
              add            IMID, %rsi
              imul           N1, %rsi
              add            Ix, %rsi
              lea            (X, %rsi, 8), %rsi
              vmovapd        (%rsi), %ymm0
              vmovapd        (%rsi, N1N2, 8), %ymm1
              vpermpd        $78, %ymm0, %ymm4
              vinsertf128    $1, %xmm1, %ymm0, %ymm0
              vinsertf128    $0, %xmm4, %ymm1, %ymm1
              vmovapd        %ymm0, (%rsi)
              vmovapd        %ymm1, (%rsi, N1N2, 8)
              vmovapd        32(%rsi, %rax, 8), %ymm0
              vmovapd        32(%rsi, %rbx, 8), %ymm1
              vpermpd        $78, %ymm0, %ymm4
              vinsertf128    $1, %xmm1, %ymm0, %ymm0
              vinsertf128    $0, %xmm4, %ymm1, %ymm1
              vmovapd        %ymm0, 32(%rsi, %rax, 8)
              vmovapd        %ymm1, 32(%rsi, %rbx, 8)
              vmovapd        32(%rsi), %ymm0
              vmovapd        32(%rsi, N1N2, 8), %ymm1
              vmovapd        (%rsi, %rax, 8), %ymm2
              vmovapd        (%rsi, %rbx, 8), %ymm3
              vpermpd        $78, %ymm0, %ymm4
              vinsertf128    $1, %xmm1, %ymm0, %ymm0
              vinsertf128    $0, %xmm4, %ymm1, %ymm1
              vpermpd        $78, %ymm2, %ymm4
              vinsertf128    $1, %xmm3, %ymm2, %ymm2
              vinsertf128    $0, %xmm4, %ymm3, %ymm3
              vmovapd        %ymm3, 32(%rsi, %rax, 8)
              vmovapd        %ymm2, 32(%rsi, %rbx, 8)
              vmovapd        32(%rsi), %ymm1
              vmovapd        32(%rsi, N1N2, 8), %ymm0
skip_diagonal:inc            IMID
              cmp            IMID, N2
              jg             main_loop
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





