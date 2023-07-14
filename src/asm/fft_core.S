#define       er0            %ymm0
#define       ei0            %ymm1
#define       er1            %ymm2
#define       ei1            %ymm3
#define       er2            %ymm4
#define       ei2            %ymm5
#define       er3            %ymm6
#define       ei3            %ymm7
#define       tr0            %ymm8
#define       tr1            %ymm9
#define       tr2            %ymm10
#define       tr3            %ymm11
#define       ti0            %ymm12
#define       ti1            %ymm13
#define       ti2            %ymm14
#define       ti3            %ymm15
#define       cos1           %ymm8
#define       sin1           %ymm9
#define       cos2           %ymm10
#define       sin2           %ymm11
#define       tr             %ymm12
#define       ti             %ymm13
#define       two            %ymm14
#define       X              %r15
#define       Wr             %r14
#define       Wi             %r13
#define       N              %r12
#define       k2             %r11
#define       level          %r10
#define       M              %r9
#define       I              %r8

              .global        fft_core2


              .text

fft_core2:    push           %rbx
              push           %r15
              push           %r13
              push           %r14
              push           %r12
              mov            %rdi, X
              mov            %rsi, N
              mov            N, %rdi
              call           get_reverse_twiddles
              mov            %rax, Wr
              mov            %rdx, Wi
              xor            level, level
              xor            k2, k2
              call           next_level
              pop            %r12
              pop            %r14
              pop            %r13
              pop            %r15
              pop            %rbx
              ret

next_level:   test           level, level
              jz             callfirst
              bsf            N, %rax
              bt             $0, %rax
              jc             call2
              cmp            $16, N
              je             call16
              call           radix4
              jmp            dispatch4
callfirst:    call           first_radix4
              jmp            dispatch4
call2:        call           radix2
              jmp            dispatch2
call16:       call           radix16
              jmp            done
dispatch2:    push           k2
              mov            N, %rdx
              shl            $3, %rdx
              shl            k2
              push           k2
              push           X
              lea            (X, %rdx), X
              inc            k2
              shr            N
              inc            level
              call           next_level
              pop            X
              pop            k2
              call           next_level
              shl            N
              dec            level
              pop            k2
              jmp            done
dispatch4:    push           k2
              mov            N, %rdx
              shl            $2, %rdx
              shl            $2, k2
              push           k2
              push           X
              lea            (X, %rdx), X
              inc            k2
              push           k2
              push           X
              lea            (X, %rdx), X
              inc            k2
              push           k2
              push           X
              lea            (X, %rdx), X
              inc            k2
              shr            $2, N
              add            $2, level
              call           next_level
              pop            X
              pop            k2
              call           next_level
              pop            X
              pop            k2
              call           next_level
              pop            X
              pop            k2
              call           next_level
              shl            $2, N
              sub            $2, level
              pop            k2
done:         ret

first_radix4: mov            N, M
              shr            $2, M
              xor            I, I
              mov            M, %rdi
              shl            $4, %rdi
              imul           $16, M, %rax
              imul           $48, M, %rbx
first_loop:   mov            I, %rdx
              shl            $4, %rdx
              lea            (X, %rdx), %rsi
              vmovapd        (%rsi), er0
              vmovapd        (%rsi, %rax), er1
              vmovapd        (%rsi, %rax, 2), er2
              vmovapd        (%rsi, %rbx), er3
              vmovapd        32(%rsi), ei0
              vmovapd        32(%rsi, %rax), ei1
              vmovapd        32(%rsi, %rax, 2), ei2
              vmovapd        32(%rsi, %rbx), ei3
              vpermpd        $177, %ymm0, %ymm8
              vpermpd        $177, %ymm1, %ymm9
              vpermpd        $177, %ymm2, %ymm10
              vpermpd        $177, %ymm3, %ymm11
              vpermpd        $177, %ymm4, %ymm12
              vpermpd        $177, %ymm5, %ymm13
              vpermpd        $177, %ymm6, %ymm14
              vpermpd        $177, %ymm7, %ymm15
              vblendpd       $5, %ymm0, %ymm9, %ymm0
              vblendpd       $10, %ymm1, %ymm8, %ymm1
              vblendpd       $5, %ymm2, %ymm11, %ymm2
              vblendpd       $10, %ymm3, %ymm10, %ymm3
              vblendpd       $5, %ymm4, %ymm13, %ymm4
              vblendpd       $10, %ymm5, %ymm12, %ymm5
              vblendpd       $5, %ymm6, %ymm15, %ymm6
              vblendpd       $10, %ymm7, %ymm14, %ymm7
              vpermpd        $216, %ymm0, %ymm0
              vpermpd        $216, %ymm1, %ymm1
              vpermpd        $216, %ymm2, %ymm2
              vpermpd        $216, %ymm3, %ymm3
              vpermpd        $216, %ymm4, %ymm4
              vpermpd        $216, %ymm5, %ymm5
              vpermpd        $216, %ymm6, %ymm6
              vpermpd        $216, %ymm7, %ymm7
              vaddpd         er2, er0, tr0
              vaddpd         ei2, ei0, ti0
              vsubpd         er2, er0, tr2
              vsubpd         ei2, ei0, ti2
              vaddpd         er3, er1, tr1
              vaddpd         ei3, ei1, ti1
              vsubpd         er3, er1, tr3
              vsubpd         ei3, ei1, ti3
              vaddpd         tr1, tr0, er0
              vaddpd         ti1, ti0, ei0
              vaddpd         ti3, tr2, er2
              vsubpd         tr3, ti2, ei2
              vsubpd         tr1, tr0, er1
              vsubpd         ti1, ti0, ei1
              vsubpd         ti3, tr2, er3
              vaddpd         tr3, ti2, ei3
              vmovapd        er0, (%rsi)
              vmovapd        er1, (%rsi, %rax)
              vmovapd        er2, (%rsi, %rax, 2)
              vmovapd        er3, (%rsi, %rbx)
              vmovapd        ei0, 32(%rsi)
              vmovapd        ei1, 32(%rsi, %rax)
              vmovapd        ei2, 32(%rsi, %rax, 2)
              vmovapd        ei3, 32(%rsi, %rbx)
              add            $4, I
              cmp            I, M
              jne            first_loop
              ret

radix4:       vmovapd        TWO, two
              mov            N, M
              shr            $2, M
              xor            I, I
              vbroadcastsd   (Wr, k2, 8), cos1
              vbroadcastsd   (Wi, k2, 8), sin1
              vmulpd         sin1, sin1, cos2
              vmulpd         cos1, sin1, sin2
              vfmsub231pd    cos1, cos1, cos2
              vfmadd231pd    sin1, cos1, sin2
              imul           $16, M, %rax
              imul           $48, M, %rbx
radix4_loop:  mov            I, %rdx
              shl            $4, %rdx
              lea            (X, %rdx), %rsi
              vmovapd        (%rsi), er0
              vmovapd        (%rsi, %rax), er1
              vmovapd        (%rsi, %rax, 2), er2
              vmovapd        (%rsi, %rbx), er3
              vmovapd        32(%rsi), ei0
              vmovapd        32(%rsi, %rax), ei1
              vmovapd        32(%rsi, %rax, 2), ei2
              vmovapd        32(%rsi, %rbx), ei3
              call           butterfly_tw4
              vmovapd        er0, (%rsi)
              vmovapd        er1, (%rsi, %rax)
              vmovapd        er3, (%rsi, %rax, 2)
              vmovapd        er2, (%rsi, %rbx)
              vmovapd        ei0, 32(%rsi)
              vmovapd        ei1, 32(%rsi, %rax)
              vmovapd        ei3, 32(%rsi, %rax, 2)
              vmovapd        ei2, 32(%rsi, %rbx)
              add            $4, I
              cmp            I, M
              jne            radix4_loop
              ret

radix2:       mov            N, M
              shr            M
              xor            I, I
              vbroadcastsd   (Wr, k2, 8), cos1
              vbroadcastsd   (Wi, k2, 8), sin1
              vmulpd         sin1, sin1, cos2
              vmulpd         cos1, sin1, sin2
              vfmsub231pd    cos1, cos1, cos2
              vfmadd231pd    sin1, cos1, sin2
              vmovapd        TWO, two
              imul           $16, M, %rax
radix2_loop:  mov            I, %rdx
              shl            $4, %rdx
              lea            (X, %rdx), %rsi
              vmovapd        (%rsi), er0
              vmovapd        (%rsi, %rax), er1
              vmovapd        32(%rsi), ei0
              vmovapd        32(%rsi, %rax), ei1
              vmovapd        er0, tr
              vmovapd        ei0, ti
              vfmadd231pd    sin2, ei1, tr
              vfnmadd231pd   sin2, er1, ti
              vfnmadd132pd   cos2, tr, er1
              vfnmadd132pd   cos2, ti, ei1
              vfmsub132pd    two, er1, er0
              vfmsub132pd    two, ei1, ei0
              vmovapd        er0, (%rsi)
              vmovapd        er1, (%rsi, %rax)
              vmovapd        ei0, 32(%rsi)
              vmovapd        ei1, 32(%rsi, %rax)
              add            $4, I
              cmp            I, M
              jne            radix2_loop
              ret


radix16:      vmovapd        (X), er0
              vmovapd        64(X), er1
              vmovapd        128(X), er2
              vmovapd        192(X), er3
              vmovapd        32(X), ei0
              vmovapd        96(X), ei1
              vmovapd        160(X), ei2
              vmovapd        224(X), ei3
              vbroadcastsd   (Wr, k2, 8), cos1
              vbroadcastsd   (Wi, k2, 8), sin1
              vmulpd         sin1, sin1, cos2
              vmulpd         cos1, sin1, sin2
              vfmsub231pd    cos1, cos1, cos2
              vfmadd231pd    sin1, cos1, sin2
              vmovapd        TWO, two
              call           butterfly_tw4
              vunpcklpd      er1, er0, %ymm8
              vunpckhpd      er1, er0, %ymm9
              vunpcklpd      er2, er3, %ymm10
              vunpckhpd      er2, er3, %ymm11
              vperm2f128     $0x20, %ymm10, %ymm8, er0
              vperm2f128     $0x20, %ymm11, %ymm9, er1
              vperm2f128     $0x31, %ymm10, %ymm8, er2
              vperm2f128     $0x31, %ymm11, %ymm9, er3
              vunpcklpd      ei1, ei0, %ymm8
              vunpckhpd      ei1, ei0, %ymm9
              vunpcklpd      ei2, ei3, %ymm10
              vunpckhpd      ei2, ei3, %ymm11
              vperm2f128     $0x20, %ymm10, %ymm8, ei0
              vperm2f128     $0x20, %ymm11, %ymm9, ei1
              vperm2f128     $0x31, %ymm10, %ymm8, ei2
              vperm2f128     $0x31, %ymm11, %ymm9, ei3
              shl            $2, k2
              vmovapd        (Wr, k2, 8), cos1
              vmovapd        (Wi, k2, 8), sin1
              vmulpd         sin1, sin1, cos2
              vmulpd         cos1, sin1, sin2
              vfmsub231pd    cos1, cos1, cos2
              vfmadd231pd    sin1, cos1, sin2
              vmovapd        TWO, two
              call           butterfly_tw4
              vmovapd        %ymm6, %ymm12
              vmovapd        %ymm4, %ymm13
              vunpcklpd      %ymm1, %ymm0, %ymm8
              vunpckhpd      %ymm1, %ymm0, %ymm9
              vunpcklpd      %ymm3, %ymm2, %ymm10
              vunpckhpd      %ymm3, %ymm2, %ymm11
              vperm2f128     $0x20, %ymm10, %ymm8, %ymm0
              vperm2f128     $0x20, %ymm11, %ymm9, %ymm2
              vperm2f128     $0x31, %ymm10, %ymm8, %ymm4
              vperm2f128     $0x31, %ymm11, %ymm9, %ymm6
              vunpcklpd      %ymm7, %ymm12, %ymm8
              vunpckhpd      %ymm7, %ymm12, %ymm9
              vunpcklpd      %ymm5, %ymm13, %ymm10
              vunpckhpd      %ymm5, %ymm13, %ymm11
              vperm2f128     $0x20, %ymm10, %ymm8, %ymm1
              vperm2f128     $0x20, %ymm11, %ymm9, %ymm3
              vperm2f128     $0x31, %ymm10, %ymm8, %ymm5
              vperm2f128     $0x31, %ymm11, %ymm9, %ymm7
              vmovapd        er0, (X)
              vmovapd        ei0, 32(X)
              vmovapd        er1, 64(X)
              vmovapd        ei1, 96(X)
              vmovapd        er2, 128(X)
              vmovapd        ei2, 160(X)
              vmovapd        er3, 192(X)
              vmovapd        ei3, 224(X)
              ret

butterfly_tw4:vmovapd        er0, tr
              vmovapd        ei0, ti
              vfmadd231pd    sin2, ei2, tr
              vfnmadd231pd   sin2, er2, ti
              vfnmadd132pd   cos2, tr, er2
              vfnmadd132pd   cos2, ti, ei2
              vfmsub132pd    two, er2, er0
              vfmsub132pd    two, ei2, ei0
              vmovapd        er1, tr
              vmovapd        ei1, ti
              vfmadd231pd    sin2, ei3, tr
              vfnmadd231pd   sin2, er3, ti
              vfnmadd132pd   cos2, tr, er3
              vfnmadd132pd   cos2, ti, ei3
              vfmsub132pd    two, er3, er1
              vfmsub132pd    two, ei3, ei1
              vmovapd        er0, tr
              vmovapd        ei0, ti
              vfmadd231pd    sin1, ei1, tr
              vfnmadd231pd   sin1, er1, ti
              vfnmadd132pd   cos1, tr, er1
              vfnmadd132pd   cos1, ti, ei1
              vfmsub132pd    two, er1, er0
              vfmsub132pd    two, ei1, ei0
              vmovapd        er2, tr
              vmovapd        ei2, ti
              vfmadd231pd    cos1, ei3, tr
              vfnmadd231pd   cos1, er3, ti
              vfmadd132pd    sin1, tr, er3
              vfmadd132pd    sin1, ti, ei3
              vfmsub132pd    two, er3, er2
              vfmsub132pd    two, ei3, ei2
              ret
