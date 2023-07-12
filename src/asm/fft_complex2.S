#define       er0            %ymm0
#define       ei0            %ymm1
#define       er1            %ymm2
#define       ei1            %ymm3
#define       er2            %ymm4
#define       ei2            %ymm5
#define       er3            %ymm6
#define       ei3            %ymm7
#define       ter0           %xmm0
#define       ter1           %xmm1
#define       ter2           %xmm2
#define       ter3           %xmm3
#define       tei0           %xmm4
#define       tei1           %xmm5
#define       tei2           %xmm6
#define       tei3           %xmm7
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
#define       half           %ymm15
#define       X              %r15
#define       Wr             %r14
#define       Wi             %r13
#define       k2rev          %r12
#define       N              %r11
#define       N2             %r10
#define       NLO            %r9
#define       nlo            %r8

              .global        fft_complex2


              .text

fft_complex2: push           %rbx
              push           %r15
              push           %r13
              push           %r14
              push           %r12
              mov            %rdi, X
              mov            %rsi, N
              mov            N, %rdi
              push           X
              push           N
              call           get_twiddles
              pop            N
              pop            X
              mov            %rax, Wr
              mov            %rdx, Wi
              mov            $1, N2
              xor            k2rev, k2rev
              call           next_level
              mov            X, %rdi
              mov            N, %rsi
              shl            %rsi
              call           scramble
              pop            %r12
              pop            %r14
              pop            %r13
              pop            %r15
              pop            %rbx
              ret

next_level:   mov            N, %rax
              shr            $2, %rax
              mov            %rax, NLO
              cmp            $16, N
              je             call16
              cmp            $32, N
              je             call32
              call           radix4
              jmp            dispatch
call16:       call           radix16
              ret
call32:       call           radix32
              ret
dispatch:     push           k2rev
              mov            NLO, %rdx
              shl            $4, %rdx
              shl            $2, k2rev
              push           k2rev
              push           X
              lea            (X, %rdx), X
              inc            k2rev
              push           k2rev
              push           X
              lea            (X, %rdx), X
              inc            k2rev
              push           k2rev
              push           X
              lea            (X, %rdx), X
              inc            k2rev
              shr            $2, N
              shl            $2, N2
              call           next_level
              pop            X
              pop            k2rev
              call           next_level
              pop            X
              pop            k2rev
              call           next_level
              pop            X
              pop            k2rev
              call           next_level
              shl            $2, N
              shr            $2, N2
              pop            k2rev
              ret

radix4:       cmp            $1, N2
              jne            not_first
              xor            nlo, nlo
              mov            NLO, %rdi
              shl            $4, %rdi
first_loop:   mov            nlo, %rdx
              shl            $4, %rdx
              lea            (X, %rdx), %rax
              lea            (%rax, %rdi), %rbx
              lea            (%rbx, %rdi), %rcx
              lea            (%rcx, %rdi), %rdx
              vmovapd        (%rax), er0
              vmovapd        (%rbx), er1
              vmovapd        (%rcx), er2
              vmovapd        (%rdx), er3
              vmovapd        32(%rax), ei0
              vmovapd        32(%rbx), ei1
              vmovapd        32(%rcx), ei2
              vmovapd        32(%rdx), ei3
              call           trans21034
              call           trans10234
              call           butterfly_04
              vmovapd        er0, (%rax)
              vmovapd        er1, (%rbx)
              vmovapd        er2, (%rcx)
              vmovapd        er3, (%rdx)
              vmovapd        ei0, 32(%rax)
              vmovapd        ei1, 32(%rbx)
              vmovapd        ei2, 32(%rcx)
              vmovapd        ei3, 32(%rdx)
              add            $4, nlo
              cmp            nlo, NLO
              jne            first_loop
              jmp            radix4_end
not_first:    xor            nlo, nlo
              vbroadcastsd   (Wr, k2rev, 8), cos1
              vbroadcastsd   (Wi, k2rev, 8), sin1
              vmulpd         sin1, sin1, cos2
              vmulpd         cos1, sin1, sin2
              vfmsub231pd    cos1, cos1, cos2
              vfmadd231pd    sin1, cos1, sin2
              imul           $16, NLO, %rax
              imul           $32, NLO, %rbx
              imul           $48, NLO, %rcx
radix4_loop:  mov            nlo, %rdx
              shl            $4, %rdx
              lea            (X, %rdx), %rsi
              vmovapd        (%rsi), er0
              vmovapd        (%rsi, %rax), er1
              vmovapd        (%rsi, %rbx), er2
              vmovapd        (%rsi, %rcx), er3
              vmovapd        32(%rsi), ei0
              vmovapd        32(%rsi, %rax), ei1
              vmovapd        32(%rsi, %rbx), ei2
              vmovapd        32(%rsi, %rcx), ei3
              call           butterfly_tw4
              vmovapd        er0, (%rsi)
              vmovapd        er1, (%rsi, %rax)
              vmovapd        er2, (%rsi, %rbx)
              vmovapd        er3, (%rsi, %rcx)
              vmovapd        ei0, 32(%rsi)
              vmovapd        ei1, 32(%rsi, %rax)
              vmovapd        ei2, 32(%rsi, %rbx)
              vmovapd        ei3, 32(%rsi, %rcx)
              add            $4, nlo
              cmp            nlo, NLO
              jne            radix4_loop
radix4_end:   ret

radix16:      vmovapd        (X), er0
              vmovapd        64(X), er1
              vmovapd        128(X), er2
              vmovapd        192(X), er3
              vmovapd        32(X), ei0
              vmovapd        96(X), ei1
              vmovapd        160(X), ei2
              vmovapd        224(X), ei3
              cmp            $1, N2
              jne            skip_trans16
              call           trans21034
              call           trans10234
skip_trans16: vbroadcastsd   (Wr, k2rev, 8), cos1
              vbroadcastsd   (Wi, k2rev, 8), sin1
              vmulpd         sin1, sin1, cos2
              vmulpd         cos1, sin1, sin2
              vfmsub231pd    cos1, cos1, cos2
              vfmadd231pd    sin1, cos1, sin2
              call           butterfly_tw4
              call           trans34201
              shl            $2, k2rev
              vmovapd        (Wr, k2rev, 8), cos1
              vmovapd        (Wi, k2rev, 8), sin1
              vmulpd         sin1, sin1, cos2
              vmulpd         cos1, sin1, sin2
              vfmsub231pd    cos1, cos1, cos2
              vfmadd231pd    sin1, cos1, sin2
              call           butterfly_tw4
              call           trans34201
              call           trans10234
              call           trans21034
              vmovapd        er0, (X)
              vmovapd        ei0, 32(X)
              vmovapd        er1, 64(X)
              vmovapd        ei1, 96(X)
              vmovapd        er2, 128(X)
              vmovapd        ei2, 160(X)
              vmovapd        er3, 192(X)
              vmovapd        ei3, 224(X)
              ret

radix32:      ret


butterfly_tw4:vmovapd        TWO, two
              vmovapd        er0, tr
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
              vmovapd        er2, tr
              vmovapd        ei2, ti
              vmovapd        er3, er2
              vmovapd        ei3, ei2
              vmovapd        tr, er3
              vmovapd        ti, ei3
              ret

butterfly_04: vaddpd         er2, er0, tr0
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
              ret

trans34201:   vunpcklpd      er1, er0, %ymm8
              vunpckhpd      er1, er0, %ymm9
              vunpcklpd      er3, er2, %ymm10
              vunpckhpd      er3, er2, %ymm11
              vperm2f128     $0x20, %ymm10, %ymm8, er0
              vperm2f128     $0x20, %ymm11, %ymm9, er1
              vperm2f128     $0x31, %ymm10, %ymm8, er2
              vperm2f128     $0x31, %ymm11, %ymm9, er3
              vunpcklpd      ei1, ei0, %ymm8
              vunpckhpd      ei1, ei0, %ymm9
              vunpcklpd      ei3, ei2, %ymm10
              vunpckhpd      ei3, ei2, %ymm11
              vperm2f128     $0x20, %ymm10, %ymm8, ei0
              vperm2f128     $0x20, %ymm11, %ymm9, ei1
              vperm2f128     $0x31, %ymm10, %ymm8, ei2
              vperm2f128     $0x31, %ymm11, %ymm9, ei3
              ret

trans10234:   vpermpd        $216, %ymm0, %ymm0
              vpermpd        $216, %ymm1, %ymm1
              vpermpd        $216, %ymm2, %ymm2
              vpermpd        $216, %ymm3, %ymm3
              vpermpd        $216, %ymm4, %ymm4
              vpermpd        $216, %ymm5, %ymm5
              vpermpd        $216, %ymm6, %ymm6
              vpermpd        $216, %ymm7, %ymm7
              ret

trans21034:   vpermpd        $177, %ymm0, %ymm8
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
              ret

              .align         32
TWO:          .double        2.0
              .double        2.0
              .double        2.0
              .double        2.0
