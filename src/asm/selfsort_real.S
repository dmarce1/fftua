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
#define       Y              %r14
#define       Wr             %r13
#define       Wi             %r12
#define       k2             %r11
#define       N2             %r10
#define       lev            %r9
#define       I              %r8
#define       STACK_SIZE     $32
#define       Wptr           -8(%rbp)
#define       N              -16(%rbp)
#define       NHI            -24(%rbp)
#define       dyn_size       -32(%rbp)

              .global        fft_selfsort_real


              .text

fft_selfsort_real:
              push           %rbp
              mov            %rsp, %rbp
              sub            STACK_SIZE, %rsp
              push           %r12
              push           %r13
              push           %r14
              push           %r15
              push           %rbx
              mov            %rdi, X
              shr            %rsi
              mov            %rsi, N
              bsf            N, %rax
              inc            %rax
              shl            $4, %rax
              sub            %rax, %rsp
              mov            %rsp, Wptr
              mov            %rax, dyn_size
              mov            %rsp, %rsi
              bsf            N, %rcx
twiddle_loop: push           %rcx
              push           %rsi
              mov            $1, %rdi
              shl            %rcx, %rdi
              call           get_forward_twiddles
              pop            %rsi
              pop            %rcx
              mov            %rcx, %rbx
              shl            %rbx
              mov            %rax, (%rsi, %rbx, 8)
              mov            %rdx, 8(%rsi, %rbx, 8)
              dec            %rcx
              test           %rcx, %rcx
              jnz            twiddle_loop
              mov            N, %rax
              lea            (X, %rax, 8), Y
              xor            lev, lev
              mov            N, N2
              shr            N2
              call           next_level
              mov            X, %rdi
              mov            N, %rsi
              shl            %rsi
              call           complex2real
              add            dyn_size, %rsp
              pop            %rbx
              pop            %r15
              pop            %r14
              pop            %r13
              pop            %r12
              mov            %rbp, %rsp
              pop            %rbp
              ret

next_level:   mov            lev, %rcx
              mov            $1, %rax
              shl            %rcx, %rax
              mov            %rax, N2
              mov            N, %rbx
              shr            $2, %rbx
              cmp            %rbx, %rax
              je             call2
              jg             return
              test           lev, lev
              jz             call0
              call           radix4_pass1
              jmp            call1
call0:        call           first_level
              jmp            call1
call2:        call           radix4_pass2
              jmp            return
call1:        mov            N, %rdx
              shr            $2, %rdx
              push           Y
              push           X
              lea            (X, %rdx, 8), X
              lea            (Y, %rdx, 8), Y
              push           Y
              push           X
              lea            (X, %rdx, 8), X
              lea            (Y, %rdx, 8), Y
              push           Y
              push           X
              lea            (X, %rdx, 8), X
              lea            (Y, %rdx, 8), Y
              add            $2, lev
              mov            N, %rax
              shr            $2, %rax
              mov            %rax, N
              call           next_level
              pop            X
              pop            Y
              call           next_level
              pop            X
              pop            Y
              call           next_level
              pop            X
              pop            Y
              call           next_level
              mov            N, %rax
              mov            %rax, N2
              shl            $2, %rax
              mov            %rax, N
              sub            $2, lev
              call           radix4_pass2
return:       ret



first_level:  mov            N, %rsi
              mov            $4, %rdi
              xor            I, I
first_loop:   lea            (X, I, 8), %rax
              lea            (%rax, %rdi, 8), %rbx
              lea            (%rbx, %rdi, 8), %rcx
              lea            (%rcx, %rdi, 8), %rdx
              vmovapd        (%rax), %ymm0
              vmovapd        (%rbx), %ymm2
              vmovapd        (%rcx), %ymm4
              vmovapd        (%rdx), %ymm6
              vmovapd        (%rax, %rsi, 8), %ymm1
              vmovapd        (%rbx, %rsi, 8), %ymm3
              vmovapd        (%rcx, %rsi, 8), %ymm5
              vmovapd        (%rdx, %rsi, 8), %ymm7
// C g0 xxx g1 g2 b3
              vpermpd        $177, %ymm0, %ymm8
              vpermpd        $177, %ymm1, %ymm9
              vpermpd        $177, %ymm2, %ymm10
              vpermpd        $177, %ymm3, %ymm11
              vpermpd        $177, %ymm4, %ymm12
              vpermpd        $177, %ymm5, %ymm13
              vpermpd        $177, %ymm6, %ymm14
              vpermpd        $177, %ymm7, %ymm15
              vblendpd       $5, %ymm0, %ymm9, %ymm0
              vblendpd       $5, %ymm2, %ymm11, %ymm2
              vblendpd       $5, %ymm4, %ymm13, %ymm4
              vblendpd       $5, %ymm6, %ymm15, %ymm6
              vblendpd       $10, %ymm1, %ymm8, %ymm1
              vblendpd       $10, %ymm3, %ymm10, %ymm3
              vblendpd       $10, %ymm5, %ymm12, %ymm5
              vblendpd       $10, %ymm7, %ymm14, %ymm7
// g1 g0 xxx C g2 b3
              vaddpd         %ymm4, %ymm0, %ymm8
              vsubpd         %ymm4, %ymm0, %ymm12
              vaddpd         %ymm5, %ymm1, %ymm9
              vsubpd         %ymm5, %ymm1, %ymm13
              vaddpd         %ymm6, %ymm2, %ymm10
              vsubpd         %ymm6, %ymm2, %ymm14
              vaddpd         %ymm7, %ymm3, %ymm11
              vsubpd         %ymm7, %ymm3, %ymm15
// g1 g0 xxx C b2 g3
              vaddpd         %ymm10, %ymm8, %ymm0
              vaddpd         %ymm11, %ymm9, %ymm4
              vaddpd         %ymm14, %ymm12, %ymm2
              vsubpd         %ymm15, %ymm13, %ymm6
              vsubpd         %ymm10, %ymm8, %ymm1
              vsubpd         %ymm11, %ymm9, %ymm5
              vsubpd         %ymm14, %ymm12, %ymm3
              vaddpd         %ymm15, %ymm13, %ymm7
              vmulpd         NONE, %ymm7, %ymm7
 // g1 g0 xxx g3 b2 C
              vunpcklpd      %ymm1, %ymm0, %ymm8
              vunpckhpd      %ymm1, %ymm0, %ymm9
              vunpcklpd      %ymm3, %ymm2, %ymm10
              vunpckhpd      %ymm3, %ymm2, %ymm11
              vperm2f128     $0x20, %ymm10, %ymm8, %ymm0
              vperm2f128     $0x20, %ymm11, %ymm9, %ymm1
              vperm2f128     $0x31, %ymm10, %ymm8, %ymm2
              vperm2f128     $0x31, %ymm11, %ymm9, %ymm3
              vunpcklpd      %ymm5, %ymm4, %ymm8
              vunpckhpd      %ymm5, %ymm4, %ymm9
              vunpcklpd      %ymm7, %ymm6, %ymm10
              vunpckhpd      %ymm7, %ymm6, %ymm11
              vperm2f128     $0x20, %ymm10, %ymm8, %ymm4
              vperm2f128     $0x20, %ymm11, %ymm9, %ymm5
              vperm2f128     $0x31, %ymm10, %ymm8, %ymm6
              vperm2f128     $0x31, %ymm11, %ymm9, %ymm7
              vmovapd        %ymm0, (%rax)
              vmovapd        %ymm2, (%rbx)
              vmovapd        %ymm4, (%rcx)
              vmovapd        %ymm6, (%rdx)
              vmovapd        %ymm1, (%rax, %rsi, 8)
              vmovapd        %ymm3, (%rbx, %rsi, 8)
              vmovapd        %ymm5, (%rcx, %rsi, 8)
              vmovapd        %ymm7, (%rdx, %rsi, 8)
// g3 g2 xxx g1 g0 C
              add            $16, I
              mov            N, %rax
              cmp            I, %rax
              jne            first_loop
              vmovapd        TWO, two
              ret

radix4_pass1: mov            %rsp, %rbx
              sub            $39, %rsp
              and            $0xffffffffffffffc0, %rsp
              mov            %rbx, (%rsp)
              mov            Wptr, %rdx
              bsf            N2, %rax
              add            $2, %rax
              shl            %rax
              mov            (%rdx, %rax, 8), Wr
              mov            8(%rdx, %rax, 8), Wi
              bsf            N, %rcx
              bsf            N2, %rax
              sub            %rax, %rcx
              sub            $4, %rcx
              mov            $1, %rdx
              shl            %rcx, %rdx
              mov            %rdx, NHI
              xor            I, I
hi_loop:      xor            k2, k2
k2_loop:      mov            $-1, %rsi
butter_loop:  inc            %rsi
// k2 + N2 * (b1 + N1 * (I + NHI * b2)
              mov            I, %rdx
              shl            $2, %rdx
              add            %rsi, %rdx
              imul           N2, %rdx
              add            k2, %rdx
              mov            N2, %rdi
              shl            $2, %rdi
              imul           NHI, %rdi
              lea            (X, %rdx, 8), %rax
              lea            (%rax, %rdi, 8), %rbx
              lea            (%rbx, %rdi, 8), %rcx
              lea            (%rcx, %rdi, 8), %rdx
              mov            Y, %rdi
              sub            X, %rdi
              vmovapd        (Wr, k2, 8), cos1
              vmovapd        (Wi, k2, 8), sin1
              vmulpd         sin1, sin1, cos2
              vmulpd         cos1, sin1, sin2
              vfmsub231pd    cos1, cos1, cos2
              vfmadd231pd    sin1, cos1, sin2
              vmovapd        (%rax), er0
              vmovapd        (%rbx), er1
              vmovapd        (%rcx), er2
              vmovapd        (%rdx), er3
              vmovapd        (%rax, %rdi), ei0
              vmovapd        (%rbx, %rdi), ei1
              vmovapd        (%rcx, %rdi), ei2
              vmovapd        (%rdx, %rdi), ei3
              call           butterfly
              cmp            $3, %rsi
              je             exit_butter
              sub            $256, %rsp
              vmovapd        %ymm0, (%rsp)
              vmovapd        %ymm1, 32(%rsp)
              vmovapd        %ymm2, 64(%rsp)
              vmovapd        %ymm3, 96(%rsp)
              vmovapd        %ymm4, 128(%rsp)
              vmovapd        %ymm5, 160(%rsp)
              vmovapd        %ymm6, 192(%rsp)
              vmovapd        %ymm7, 224(%rsp)
              jmp            butter_loop
exit_butter:  inc            %rsi
              mov            N2, %rdi
// k2 + N2 * (b1 + N1 * (I + NHI * b2)
store_loop:   dec            %rsi
              mov            %rsi, %rdx
              imul           NHI, %rdx
              add            I, %rdx
              shl            $2, %rdx
              imul           N2, %rdx
              add            k2, %rdx
              mov            N2, %rdi
              lea            (X, %rdx, 8), %rax
              lea            (%rax, %rdi, 8), %rbx
              lea            (%rbx, %rdi, 8), %rcx
              lea            (%rcx, %rdi, 8), %rdx
              mov            Y, %rdi
              sub            X, %rdi
              shr            $3, %rdi
              vmovapd        er0, (%rax)
              vmovapd        er1, (%rbx)
              vmovapd        er3, (%rcx)
              vmovapd        er2, (%rdx)
              vmovapd        ei0, (%rax, %rdi, 8)
              vmovapd        ei1, (%rbx, %rdi, 8)
              vmovapd        ei3, (%rcx, %rdi, 8)
              vmovapd        ei2, (%rdx, %rdi, 8)
              lea            (%rcx, %rbx, 2), %rcx
              sub            X, %rcx
              shr            $3, %rcx
              test           %rsi, %rsi
              jz             exit_store
              vmovapd        (%rsp), %ymm0
              vmovapd        32(%rsp), %ymm1
              vmovapd        64(%rsp), %ymm2
              vmovapd        96(%rsp), %ymm3
              vmovapd        128(%rsp), %ymm4
              vmovapd        160(%rsp), %ymm5
              vmovapd        192(%rsp), %ymm6
              vmovapd        224(%rsp), %ymm7
              add            $256, %rsp
              jmp            store_loop
exit_store:   add            $4, k2
              cmp            k2, N2
              jne            k2_loop
              inc            I
              cmp            I, NHI
              jne            hi_loop
              mov            (%rsp), %rsp
              ret

radix4_pass2: bsf            N, %rax
              shl            %rax
              mov            Wptr, %rdx
              mov            (%rdx, %rax, 8), Wr
              mov            8(%rdx, %rax, 8), Wi
              xor            k2, k2
              mov            Y, %rdi
              sub            X, %rdi
              shr            $3, %rdi
radix4_loop2: vmovapd        (Wr, k2, 8), cos1
              vmovapd        (Wi, k2, 8), sin1
              vmulpd         sin1, sin1, cos2
              vmulpd         cos1, sin1, sin2
              vfmsub231pd    cos1, cos1, cos2
              vfmadd231pd    sin1, cos1, sin2
              lea            (X, k2, 8), %rax
              lea            (%rax, N2, 8), %rbx
              lea            (%rbx, N2, 8), %rcx
              lea            (%rcx, N2, 8), %rdx
              vmovapd        (%rax), er0
              vmovapd        (%rbx), er1
              vmovapd        (%rcx), er2
              vmovapd        (%rdx), er3
              vmovapd        (%rax, %rdi, 8), ei0
              vmovapd        (%rbx, %rdi, 8), ei1
              vmovapd        (%rcx, %rdi, 8), ei2
              vmovapd        (%rdx, %rdi, 8), ei3
              call           butterfly
              vmovapd        er0, (%rax)
              vmovapd        er1, (%rbx)
              vmovapd        er3, (%rcx)
              vmovapd        er2, (%rdx)
              vmovapd        ei0, (%rax, %rdi, 8)
              vmovapd        ei1, (%rbx, %rdi, 8)
              vmovapd        ei3, (%rcx, %rdi, 8)
              vmovapd        ei2, (%rdx, %rdi, 8)
              add            $4, k2
              cmp            k2, N2
              jne            radix4_loop2
              ret

butterfly:    vmovapd        er0, tr
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
