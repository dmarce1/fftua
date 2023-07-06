#define       X              %r15
#define       N              %r14
#define       NMID           %r13
#define       ihi            %r12
#define       imid           %r11
#define       ilo            %r10
#define       logSIMD_SIZE   $2
#define       SIMD_SIZE      $4



              .global        transpose_out


transpose_out:push           %rbx
              push           %r15
              push           %r14
              push           %r13
              push           %r12
              mov            %rdi, X
              mov            %rsi, N
              mov            N, NMID
              shr            logSIMD_SIZE, NMID
              shr            logSIMD_SIZE, NMID
              mov            NMID, %rax
              shl            logSIMD_SIZE, %rax
              imul           $2, %rax, %rbx
              imul           $3, %rax, %rcx
              xor            imid, imid
L3000:        mov            imid, %rsi
              shl            logSIMD_SIZE, %rsi
              lea            (X, %rsi, 8), %rsi
              vmovapd        (%rsi), %ymm0
              vmovapd        (%rsi, %rax, 8), %ymm1
              vmovapd        (%rsi, %rbx, 8), %ymm2
              vmovapd        (%rsi, %rcx, 8), %ymm3
              vshufpd        $0, %ymm1, %ymm0, %ymm4
              vshufpd        $15, %ymm1, %ymm0, %ymm5
              vshufpd        $0, %ymm3, %ymm2, %ymm6
              vshufpd        $15, %ymm3, %ymm2, %ymm7
              vinsertf128    $1, %xmm6, %ymm4, %ymm0
              vinsertf128    $1, %xmm7, %ymm5, %ymm1
              vpermpd        $78, %ymm4, %ymm4
              vpermpd        $78, %ymm5, %ymm5
              vinsertf128    $0, %xmm4, %ymm6, %ymm2
              vinsertf128    $0, %xmm5, %ymm7, %ymm3
              vmovapd        %ymm0, (%rsi)
              vmovapd        %ymm1, (%rsi, %rax, 8)
              vmovapd        %ymm2, (%rsi, %rbx, 8)
              vmovapd        %ymm3, (%rsi, %rcx, 8)
              inc            imid
              cmp            imid, NMID
              jne            L3000
              pop            %r12
              pop            %r13
              pop            %r14
              pop            %r15
              pop            %rbx
              ret
