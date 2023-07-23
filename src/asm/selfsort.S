#define       X              %r15
#define       N              %r14
#define       Wr             %r13
#define       Wi             %r12
#define       M              %r11
#define       I              %r10
#define       Y              %r9
#define       k2             %r8
#define       er0            %ymm0
#define       ei0            %ymm1
#define       er1            %ymm2
#define       ei1            %ymm3
#define       er2            %ymm4
#define       ei2            %ymm5
#define       er3            %ymm6
#define       ei3            %ymm7
#define       tr             %ymm8
#define       ti             %ymm9
#define       two            %ymm10
#define       cos1           %ymm11
#define       sin1           %ymm12
#define       cos2           %ymm13
#define       sin2           %ymm14
#define       STACK_SIZE     $40
#define       Wptrs          -8(%rbp)
#define       Wptrs_sz       -16(%rbp)
#define       N2             -24(%rbp)
#define       level          -32(%rbp)
#define       J              -40(%rbp)


              .global        fft_selfsort


              .text

fft_selfsort: push           %rbp
              mov            %rsp, %rbp
              sub            STACK_SIZE, %rsp
              push           %rbx
              push           %r15
              push           %r14
              push           %r13
              push           %r12
              mov            %rdi, X
              mov            %rsi, N
              bsr            N, %rax
              inc            %rax
              shl            $4, %rax
              sub            %rax, %rsp
              mov            %rax, Wptrs_sz
              mov            %rsp, Wptrs
              mov            $1, %rcx
twiddle_loop: push           %rcx
              mov            $1, %rdi
              shl            %rcx, %rdi
              call           get_forward_twiddles
              pop            %rcx
              mov            %rcx, %rdi
              shl            %rdi
              mov            Wptrs, %rsi
              mov            %rax, (%rsi, %rdi, 8)
              mov            %rdx, 8(%rsi, %rdi, 8)
              inc            %rcx
              bsr            N, %rax
              cmp            %rcx, %rax
              jge            twiddle_loop
              lea            (X, N, 8), Y
              xor            %rax, %rax
              mov            %rax, level
              call           pass0
              call           pass1
              add            Wptrs_sz, %rsp
              pop            %r12
              pop            %r13
              pop            %r14
              pop            %r15
              pop            %rbx
              mov            %rbp, %rsp
              pop            %rbp
              ret

pass0:        cmp            $0, level
              je             call0
              bsr            N, %rax
              sub            level, %rax
              sub            $3, %rax
              cmp            $0, %rax
              jle            call1
              call           butterflies4_and_transpose
              jmp            continue
call0:        call           first_butterflies4
              push           X
              push           Y
              lea            (X, N, 4), X
              lea            (Y, N, 4), Y
              shr            N
              mov            level, %rax
              inc            %rax
              mov            %rax, level
              call           pass0
              pop            Y
              pop            X
              call           pass0
              shl            N
              mov            level, %rax
              dec            %rax
              mov            %rax, level
call1:        ret
continue:     push           X
              push           Y
              lea            (X, N, 2), X
              lea            (Y, N, 2), Y
              push           X
              push           Y
              lea            (X, N, 2), X
              lea            (Y, N, 2), Y
              push           X
              push           Y
              lea            (X, N, 2), X
              lea            (Y, N, 2), Y
              shr            $2, N
              mov            level, %rax
              add            $2, %rax
              mov            %rax, level
              call           pass0
              pop            Y
              pop            X
              call           pass0
              pop            Y
              pop            X
              call           pass0
              pop            Y
              pop            X
              call           pass0
              shl            $2, N
              mov            level, %rax
              sub            $2, %rax
              mov            %rax, level
skip2:        ret

pass1:        bsr            N, %rax
              sub            level, %rax
              dec            %rax
              cmp            $0, %rax
              jle            done1
              cmp            $0, level
              jne            fours
              push           X
              push           Y
              lea            (X, N, 4), X
              lea            (Y, N, 4), Y
              shr            N
              mov            level, %rax
              inc            %rax
              mov            %rax, level
              call           pass1
              pop            Y
              pop            X
              call           pass1
              shl            N
              mov            level, %rax
              dec            %rax
              mov            %rax, level
              call           butterflies2
              ret
fours:        push           X
              push           Y
              lea            (X, N, 2), X
              lea            (Y, N, 2), Y
              push           X
              push           Y
              lea            (X, N, 2), X
              lea            (Y, N, 2), Y
              push           X
              push           Y
              lea            (X, N, 2), X
              lea            (Y, N, 2), Y
              shr            $2, N
              mov            level, %rax
              add            $2, %rax
              mov            %rax, level
              call           pass1
              pop            Y
              pop            X
              call           pass1
              pop            Y
              pop            X
              call           pass1
              pop            Y
              pop            X
              call           pass1
              shl            $2, N
              mov            level, %rax
              sub            $2, %rax
              mov            %rax, level
              call           butterflies4
done1:        ret

first_butterflies4:
              mov            N, M
              shr            M
              xor            I, I
first_loop:   lea            (X, I, 8), %rax
              lea            (%rax, M, 8), %rbx
              lea            (%rbx, M, 8), %rcx
              lea            (%rcx, M, 8), %rdx
              vmovapd        (%rax), %ymm0
              vmovapd        (%rbx), %ymm1
              vmovapd        (%rcx), %ymm2
              vmovapd        (%rdx), %ymm3
              vaddpd         %ymm2, %ymm0, %ymm4
              vsubpd         %ymm2, %ymm0, %ymm6
              vaddpd         %ymm3, %ymm1, %ymm5
              vsubpd         %ymm3, %ymm1, %ymm7
              vpermpd        $177, %ymm7, %ymm7
              vmulpd         ODD_NEG, %ymm7, %ymm7
              vaddpd         %ymm5, %ymm4, %ymm0
              vsubpd         %ymm5, %ymm4, %ymm2
              vaddpd         %ymm7, %ymm6, %ymm1
              vsubpd         %ymm7, %ymm6, %ymm3
              vunpcklpd      %ymm1, %ymm0, %ymm8
              vunpckhpd      %ymm1, %ymm0, %ymm9
              vunpcklpd      %ymm3, %ymm2, %ymm10
              vunpckhpd      %ymm3, %ymm2, %ymm11
              vperm2f128     $0x20, %ymm10, %ymm8, %ymm0
              vperm2f128     $0x20, %ymm11, %ymm9, %ymm2
              vperm2f128     $0x31, %ymm10, %ymm8, %ymm1
              vperm2f128     $0x31, %ymm11, %ymm9, %ymm3
              vmovapd        %ymm0, (%rax)
              vmovapd        %ymm1, (%rbx)
              vmovapd        %ymm2, (%rcx)
              vmovapd        %ymm3, (%rdx)
              add            $4, I
              cmp            I, M
              jne            first_loop
              ret

butterflies4_and_transpose:
              mov            level, %rax
              add            $3, %rax
              shl            %rax
              mov            Wptrs, %rdx
              mov            (%rdx, %rax, 8), Wr
              mov            8(%rdx, %rax, 8), Wi
              mov            %rsp, %rax
              and            $0xffffffffffffffc0, %rsp
              push           %rax
              push           %rax
              push           %rax
              push           %rax
              mov            level, %rcx
              inc            %rcx
              mov            $1, %rax
              shl            %rcx, %rax
              mov            %rax, N2
              bsr            N, %rax
              sub            %rax, %rcx
              neg            %rcx
              sub            $4, %rcx
              mov            $1, M
              shl            %rcx, M
              xor            I, I
outer_loop:   xor            k2, k2
inner_loop:   vmovapd        (Wr, k2, 8), cos1
              vmovapd        (Wi, k2, 8), sin1
              vmulpd         sin1, sin1, cos2
              vmulpd         cos1, sin1, sin2
              vfmsub231pd    cos1, cos1, cos2
              vfmadd231pd    sin1, cos1, sin2
              xor            %rax, %rax
              mov            N, %rdi
              shl            %rdi
butter_loop:  mov            I, %rsi
              shl            $2, %rsi
              add            %rax, %rsi
              imul           N2, %rsi
              add            k2, %rsi
              push           %rax
              lea            (X, %rsi, 8), %rax
              lea            (%rax, %rdi, 2), %rbx
              lea            (Y, %rsi, 8), %rcx
              lea            (%rcx, %rdi, 2), %rdx
              vmovapd        (%rax), er0
              vmovapd        (%rax, %rdi), er1
              vmovapd        (%rbx), er2
              vmovapd        (%rbx, %rdi), er3
              vmovapd        (%rcx), ei0
              vmovapd        (%rcx, %rdi), ei1
              vmovapd        (%rdx), ei2
              vmovapd        (%rdx, %rdi), ei3
              pop            %rax
              call           butterfly4
              inc            %rax
              cmp            $4, %rax
              je             exit_butter
              sub            $256, %rsp
              vmovapd        er0, (%rsp)
              vmovapd        er1, 32(%rsp)
              vmovapd        er2, 64(%rsp)
              vmovapd        er3, 96(%rsp)
              vmovapd        ei0, 128(%rsp)
              vmovapd        ei1, 160(%rsp)
              vmovapd        ei2, 192(%rsp)
              vmovapd        ei3, 224(%rsp)
              jmp            butter_loop
exit_butter:  mov            $3, %rax
              mov            N2, %rdi
              shl            $3, %rdi
store_loop:   mov            %rax, %rsi
              imul           M, %rsi
              add            I, %rsi
              shl            $2, %rsi
              imul           N2, %rsi
              add            k2, %rsi
              push           %rax
              lea            (X, %rsi, 8), %rax
              lea            (%rax, %rdi, 2), %rbx
              lea            (Y, %rsi, 8), %rcx
              lea            (%rcx, %rdi, 2), %rdx
              vmovapd        er0, (%rax)
              vmovapd        er3, (%rax, %rdi)
              vmovapd        er1, (%rbx)
              vmovapd        er2, (%rbx, %rdi)
              vmovapd        ei0, (%rcx)
              vmovapd        ei3, (%rcx, %rdi)
              vmovapd        ei1, (%rdx)
              vmovapd        ei2, (%rdx, %rdi)
              pop            %rax
              dec            %rax
              cmp            $-1, %rax
              je             exit_store
              vmovapd        (%rsp), er0
              vmovapd        32(%rsp), er1
              vmovapd        64(%rsp), er2
              vmovapd        96(%rsp), er3
              vmovapd        128(%rsp), ei0
              vmovapd        160(%rsp), ei1
              vmovapd        192(%rsp), ei2
              vmovapd        224(%rsp), ei3
              add            $256, %rsp
              jmp            store_loop
exit_store:   add            $4, k2
              cmp            k2, N2
              jne            inner_loop
              inc            I
              cmp            I, M
              jne            outer_loop
              pop            %rax
              pop            %rax
              pop            %rax
              pop            %rax
              mov            %rax, %rsp
              ret

butterflies4: bsr            N, %rax
              shl            %rax
              mov            Wptrs, %rdx
              mov            (%rdx, %rax, 8), Wr
              mov            8(%rdx, %rax, 8), Wi
              mov            N, %rax
              shr            $2, %rax
              mov            %rax, N2
              xor            k2, k2
butter_loop2: vmovapd        (Wr, k2, 8), cos1
              vmovapd        (Wi, k2, 8), sin1
              vmulpd         sin1, sin1, cos2
              vmulpd         cos1, sin1, sin2
              vfmsub231pd    cos1, cos1, cos2
              vfmadd231pd    sin1, cos1, sin2
              lea            (X, k2, 8), %rax
              lea            (%rax, N, 4), %rbx
              lea            (Y, k2, 8), %rcx
              lea            (%rcx, N, 4), %rdx
              vmovapd        (%rax), er0
              vmovapd        (%rbx), er2
              vmovapd        (%rax, N, 2), er1
              vmovapd        (%rbx, N, 2), er3
              vmovapd        (%rcx), ei0
              vmovapd        (%rdx), ei2
              vmovapd        (%rcx, N, 2), ei1
              vmovapd        (%rdx, N, 2), ei3
              call           butterfly4
              vmovapd        er0, (%rax)
              vmovapd        er3, (%rax, N, 2)
              vmovapd        er1, (%rbx)
              vmovapd        er2, (%rbx, N, 2)
              vmovapd        ei0, (%rcx)
              vmovapd        ei3, (%rcx, N, 2)
              vmovapd        ei1, (%rdx)
              vmovapd        ei2, (%rdx, N, 2)
              add            $4, k2
              cmp            k2, N2
              jne            butter_loop2
              ret

butterflies2: bsr            N, %rax
              shl            %rax
              mov            Wptrs, %rdx
              mov            (%rdx, %rax, 8), Wr
              mov            8(%rdx, %rax, 8), Wi
              mov            N, %rax
              shr            %rax
              mov            %rax, N2
              xor            k2, k2
butter_loop3: vmovapd        (Wr, k2, 8), cos1
              vmovapd        (Wi, k2, 8), sin1
              lea            (X, k2, 8), %rax
              lea            (Y, k2, 8), %rcx
              lea            (%rax, N, 4), %rbx
              lea            (%rcx, N, 4), %rdx
              vmovapd        (%rax), er0
              vmovapd        (%rbx), er1
              vmovapd        (%rcx), ei0
              vmovapd        (%rdx), ei1
              vmulpd         sin1, ei1, tr
              vmulpd         cos1, ei1, ti
              vfmsub231pd    cos1, er1, tr
              vfmadd231pd    sin1, er1, ti
              vaddpd         tr, er0, er2
              vsubpd         tr, er0, er3
              vaddpd         ti, ei0, ei2
              vsubpd         ti, ei0, ei3
              vmovapd        er2, (%rax)
              vmovapd        er3, (%rbx)
              vmovapd        ei2, (%rcx)
              vmovapd        ei3, (%rdx)
              add            $4, k2
              cmp            k2, N2
              jne            butter_loop3
              ret

butterfly4:   vmovapd        TWO, two
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
              ret

              .align         32
TWO:          .double        2.0
              .double        2.0
              .double        2.0
              .double        2.0
ODD_NEG:      .double        1.0
              .double        -1.0
              .double        1.0
              .double        -1.0
