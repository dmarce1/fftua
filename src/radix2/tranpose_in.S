#define       X              %r8
#define       M              %r9
#define       Ix             %r10
#define       Iy             %r11
#define       imid           %r12
#define       ihi            %r13
#define       N1             %r14
#define       D              %r15
#define       M0             $8
#define       STACK_SIZE     $16
#define       N              -8(%rbp)
#define       N2             -16(%rbp)

              .global        transpose_in


              .text


transpose_in: push           %rbp
              mov            %rsp, %rbp
              sub            STACK_SIZE, %rsp
              push           %r12
              push           %r13
              push           %r14
              push           %r15
              push           %rbx
              mov            %rdi, X
              mov            %rsi, N
              bsf            N, %rcx
              shr            %rcx
              mov            $1, N1
              sub            $2, %rcx
              shl            %rcx, N1
              bsf            N, %rax
              and            $1, %rax
              inc            %rax
              mov            %rax, N2
              mov            N1, %rax
              imul           N2, %rax
              mov            %rax, D
              xor            Ix, Ix
              xor            Iy, Iy
              mov            N1, M
              call           recurse
              pop            %rbx
              pop            %r15
              pop            %r14
              pop            %r13
              pop            %r12
              mov            %rbp, %rsp
              pop            %rbp
              ret
recurse:      cmp            Ix, Iy
              jl             done
              cmp            M0, M
              jg             further
              cmp            Ix, Iy
              je             eq
              xor            ihi, ihi
loop_hi:      xor            imid, imid
loop_mid:     mov            ihi, %rsi
              imul           N1, %rsi
              mov            %rsi, %rdi
              add            Ix, %rsi
              add            Iy, %rdi
              imul           N2, %rsi
              imul           N2, %rdi
              add            imid, %rsi
              add            imid, %rdi
              imul           N1, %rsi
              imul           N1, %rdi
              add            Iy, %rsi
              add            Ix, %rdi
              shl            $4, %rsi
              shl            $4, %rdi
              lea            (X, %rsi, 8), %rsi
              lea            (X, %rdi, 8), %rdi
              xor            %rax, %rax
L100:         xor            %rbx, %rbx
L200:         mov            %rax, %rcx
              mov            %rbx, %rdx
              imul           D, %rcx
              imul           D, %rdx
              add            %rbx, %rcx
              add            %rax, %rdx
              shl            $2, %rcx
              shl            $2, %rdx
              vmovapd        (%rsi, %rcx, 8), %ymm0
              vmovapd        (%rdi, %rdx, 8), %ymm1
              vmovapd        %ymm1, (%rsi, %rcx, 8)
              vmovapd        %ymm0, (%rdi, %rdx, 8)
              inc            %rbx
              cmp            M0, %rbx
              jl             L200
              inc            %rax
              cmp            M0, %rax
              jne            L100
              inc            imid
              cmp            imid, N2
              jg             loop_mid
              inc            ihi
              cmp            $4, ihi
              jg             loop_hi
eq:           xor            ihi, ihi
loop_hi_eq:   xor            imid, imid
loop_mid_eq:  mov            ihi, %rsi
              imul           N1, %rsi
              add            Ix, %rsi
              imul           N2, %rsi
              add            imid, %rsi
              imul           N1, %rsi
              add            Iy, %rsi
              shl            $4, %rsi
              lea            (X, %rsi, 8), %rsi
              xor            %rax, %rax
L300:         mov            %rax, %rbx
              inc            %rbx
L400:         mov            %rax, %rcx
              mov            %rbx, %rdx
              imul           D, %rcx
              imul           D, %rdx
              add            %rbx, %rcx
              add            %rax, %rdx
              shl            $2, %rcx
              shl            $2, %rdx
              vmovapd        (%rsi, %rcx, 8), %ymm0
              vmovapd        (%rsi, %rdx, 8), %ymm1
              vmovapd        %ymm1, (%rsi, %rcx, 8)
              vmovapd        %ymm0, (%rsi, %rdx, 8)
              inc            %rbx
              cmp            M0, %rbx
              jl             L400
              inc            %rax
              cmp            M0, %rax
              jne            L300
              inc            imid
              cmp            imid, N2
              jg             loop_mid_eq
              inc            ihi
              cmp            $4, ihi
              jg             loop_hi_eq
              jmp            done
further:      shr            M
              call           recurse
              add            M, Iy
              call           recurse
              sub            M, Iy
              add            M, Ix
              call           recurse
              add            M, Iy
              call           recurse
              sub            M, Iy
              sub            M, Ix
              shl            M
done:         ret





