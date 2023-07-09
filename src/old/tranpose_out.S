#define       Xin            %r15
#define       Xout           %r14
#define       N1             %r13
#define       N2             %r12
#define       loc            %r11
#define       n2             %r10

              .global        transpose_out

              .text

transpose_out:push           %r12
              push           %r13
              push           %r14
              push           %r15
              push           %rbx
              mov            %rdi, Xin
              mov            %rsi, Xout
              mov            %rdx, N1
              mov            %rcx, N2
              mov            %r8, loc
              imul           $2, N2, %rax
              imul           $3, N2, %rbx
              imul           $2, N1, %rcx
              imul           $3, N1, %rdx
              xor            n2, n2
main_loop:    mov            loc, %rdi
              imul           N1, %rdi
              add            n2, %rdi
              mov            n2, %rsi
              lea            (Xin, %rsi, 8), %rsi
              lea            (Xout, %rdi, 8), %rdi
              vmovapd        (%rsi), %ymm0
              vmovapd        (%rsi, N2, 8), %ymm1
              vmovapd        (%rsi, %rax, 8), %ymm2
              vmovapd        (%rsi, %rbx, 8), %ymm3
              vunpcklpd      %ymm1, %ymm0, %ymm4
              vunpckhpd      %ymm1, %ymm0, %ymm5
              vunpcklpd      %ymm3, %ymm2, %ymm6
              vunpckhpd      %ymm3, %ymm2, %ymm7
              vperm2f128     $0x20, %ymm6, %ymm4, %ymm0
              vperm2f128     $0x20, %ymm7, %ymm5, %ymm1
              vperm2f128     $0x31, %ymm6, %ymm4, %ymm2
              vperm2f128     $0x31, %ymm7, %ymm5, %ymm3
              vmovapd        %ymm0, (%rdi)
              vmovapd        %ymm1, (%rdi, N1, 8)
              vmovapd        %ymm2, (%rdi, %rcx, 8)
              vmovapd        %ymm3, (%rdi, %rdx, 8)
              add            $4, n2
              cmp            n2, N2
              jne            main_loop
              pop            %rbx
              pop            %r15
              pop            %r14
              pop            %r13
              pop            %r12
              ret
