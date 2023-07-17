#define       er0            %ymm0
#define       ei0            %ymm1
#define       er1            %ymm2
#define       ei1            %ymm3
#define       er2            %ymm4
#define       ei2            %ymm5
#define       er3            %ymm6
#define       ei3            %ymm7
#define       cos1           %ymm8
#define       sin1           %ymm9
#define       cos2           %ymm10
#define       sin2           %ymm11
#define       tr             %ymm12
#define       ti             %ymm13
#define       ter0           %xmm0
#define       tei0           %xmm1
#define       ter1           %xmm2
#define       tei1           %xmm3
#define       ter2           %xmm4
#define       tei2           %xmm5
#define       ter3           %xmm6
#define       tei3           %xmm7
#define       tcos1          %xmm8
#define       tsin1          %xmm9
#define       tcos2          %xmm10
#define       tsin2          %xmm11
#define       ttr            %xmm12
#define       tti            %xmm13
#define       X              %r15
#define       Wr             %r14
#define       Wi             %r13
#define       k2             %r12
#define       N2             %r11
#define       NHI            %r10
#define       I              %r9
#define       N              %r8
#define       STACK_SIZE     $32
#define       Wptr           -8(%rbp)
#define       lev            -16(%rbp)
#define       dyn_size       -24(%rbp)

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
              push           %r8
              call           get_forward_twiddles
              pop            %r8
              pop            %rsi
              pop            %rcx
              mov            %rcx, %rbx
              shl            %rbx
              mov            %rax, (%rsi, %rbx, 8)
              mov            %rdx, 8(%rsi, %rbx, 8)
              dec            %rcx
              test           %rcx, %rcx
              jnz            twiddle_loop
              xor            %rax, %rax
              mov            %rax, lev
              mov            N, N2
              mov            $1, N2
              call           next_level
              add            dyn_size, %rsp
              pop            %rbx
              pop            %r15
              pop            %r14
              pop            %r13
              pop            %r12
              mov            %rbp, %rsp
              pop            %rbp
              ret


next_level:   bsf            N, %rax
              mov            lev, %rcx
              test           %rcx, %rcx
              jz             call_first
              sub            %rcx, %rax
              cmp            $0, %rcx
              jle            call0
              cmp            $2, %rcx
              jl             call1
              je             call2
              cmp            $3, %rcx
              je             call3
              call           radix4_selfsort
              jmp            continue
call_first:   call           radix4_first
              jmp            continue
call0:        ret
call1:        hlt
call2:        mov            N, N2
              shr            $2, N2
              call           radix4
              ret
call3:        hlt
continue:     mov            N, %rdx
              shr            $2, %rdx
              push           X
              lea            (X, %rdx, 8), X
              push           X
              lea            (X, %rdx, 8), X
              push           X
              lea            (X, %rdx, 8), X
              add            $2, lev
              shr            $2, N
              call           next_level
              pop            X
              call           next_level
              pop            X
              call           next_level
              pop            X
              call           next_level
              mov            N, N2
              shl            $2, N
              sub            $2, lev
              call           radix4
              ret

radix4_first: mov            N, %rdi
			  shr            $2, %rdi
			  xor            I, I
first_loop:   lea            (X, I, 8), %rax
              lea            (%rax, %rdi, 8), %rbx
              lea            (%rbx, %rdi, 8), %rcx
              lea            (%rcx, %rdi, 8), %rdx
              vmovapd        (%rax), er0
              vmovapd        (%rbx), er1
              vmovapd        (%rcx), er2
              vmovapd        (%rdx), er3
              call           butter0_4
              vunpcklpd      er1, er0, %ymm8
              vunpckhpd      er1, er0, %ymm9
              vunpcklpd      er3, er2, %ymm10
              vunpckhpd      er3, er2, %ymm11
              vperm2f128     $0x20, %ymm10, %ymm8, er0
              vperm2f128     $0x20, %ymm11, %ymm9, er1
              vperm2f128     $0x31, %ymm10, %ymm8, er2
              vperm2f128     $0x31, %ymm11, %ymm9, er3
              vmovapd        er0, (%rax)
              vmovapd        er1, (%rbx)
              vmovapd        er2, (%rcx)
              vmovapd        er3, (%rdx)
              add            $4, I
              cmp            I, %rdi
              jne            first_loop
              ret

radix4_selfsort:
              mov            %rsp, %rbx
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
              mov            I, %rdx
              shl            $2, %rdx
              add            %rsi, %rdx
              imul           N2, %rdx
              add            k2, %rdx
              mov            N2, %rdi
              shl            $2, %rdi
              imul           NHI, %rdi
              lea            (X, %rdx, 8), %rax


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
store_loop:   dec            %rsi
              mov            %rsi, %rdx
              imul           NHI, %rdx
              add            I, %rdx
              shl            $2, %rdx
              imul           N2, %rdx
              add            k2, %rdx
              mov            N2, %rdi
              lea            (X, %rdx, 8), %rax



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

radix4:       bsf            N, %rax
              shl            %rax
              mov            Wptr, %rdx
              mov            (%rdx, %rax, 8), Wr
              mov            8(%rdx, %rax, 8), Wi
              lea            (X), %rax
              lea            (%rax, N2, 8), %rbx
              lea            (%rbx, N2, 8), %rcx
              lea            (%rcx, N2, 8), %rdx
              vmovq          (%rax), ter0
              vmovq          (%rbx), ter1
              vmovq          (%rcx), ter2
              vmovq          (%rdx), ter3
              call           butter0_1
              vmovq          ter0, (%rax)
              vmovq          ter1, (%rbx)
              vmovq          ter2, (%rcx)
              vmovq          ter3, (%rdx)
              cmp            $2, N2
              jl             radix4_end
              lea            (X, N2, 4), %rax
              lea            (%rax, N2, 8), %rbx
              lea            (%rbx, N2, 8), %rcx
              lea            (%rcx, N2, 8), %rdx
              vmovq          (%rax), ter0
              vmovq          (%rbx), ter1
              vmovq          (%rcx), ter2
              vmovq          (%rdx), ter3
              call           butterNy_1
              vmovq          ter0, (%rax)
              vmovq          ter1, (%rbx)
              vmovq          ter2, (%rcx)
              vmovq          ter3, (%rdx)
              cmp            $4, N2
              jl             radix4_end
              vmovq          8(Wr), tcos1
              vmovq          8(Wi), tsin1
              vmulsd         tsin1, tsin1, tcos2
              vmulsd         tcos1, tsin1, tsin2
              vfmsub231sd    tcos1, tcos1, tcos2
              vfmadd231sd    tsin1, tcos1, tsin2
              mov            N2, %rdi
              sub            $2, %rdi
              lea            8(X), %rax
              lea            (%rax, N2, 8), %rbx
              lea            (%rbx, N2, 8), %rcx
              lea            (%rcx, N2, 8), %rdx
              vmovq          (%rax), ter0
              vmovq          (%rbx), ter1
              vmovq          (%rcx), ter2
              vmovq          (%rdx), ter3
              vmovq          (%rax, %rdi, 8), tei0
              vmovq          (%rbx, %rdi, 8), tei1
              vmovq          (%rcx, %rdi, 8), tei2
              vmovq          (%rdx, %rdi, 8), tei3
              call           butter_1
              vmovq          ter0, (%rax)
              vmovq          ter2, (%rax, %rdi, 8)
              vmovq          ter3, (%rbx)
              vmovq          ter1, (%rbx, %rdi, 8)
              vmovq          tei1, (%rcx)
              vmovq          tei3, (%rcx, %rdi, 8)
              vmovq          tei2, (%rdx)
              vmovq          tei0, (%rdx, %rdi, 8)
              cmp            $8, N2
              jl             radix4_end
              vmovapd        16(Wr), tcos1
              vmovapd        16(Wi), tsin1
              vmulpd         tsin1, tsin1, tcos2
              vmulpd         tcos1, tsin1, tsin2
              vfmsub231pd    tcos1, tcos1, tcos2
              vfmadd231pd    tsin1, tcos1, tsin2
              mov            N2, %rdi
              sub            $5, %rdi
              lea            16(X), %rax
              lea            (%rax, N2, 8), %rbx
              lea            (%rbx, N2, 8), %rcx
              lea            (%rcx, N2, 8), %rdx
              vmovapd        (%rax), ter0
              vmovapd        (%rbx), ter1
              vmovapd        (%rcx), ter2
              vmovapd        (%rdx), ter3
              vmovupd        (%rax, %rdi, 8), tei0
              vmovupd        (%rbx, %rdi, 8), tei1
              vmovupd        (%rcx, %rdi, 8), tei2
              vmovupd        (%rdx, %rdi, 8), tei3
              vpshufb        PERMUTE, tei0, tei0
              vpshufb        PERMUTE, tei1, tei1
              vpshufb        PERMUTE, tei2, tei2
              vpshufb        PERMUTE, tei3, tei3
              call           butter_2
              vpshufb        PERMUTE, ter2, ter2
              vpshufb        PERMUTE, ter1, ter1
              vpshufb        PERMUTE, tei3, tei3
              vpshufb        PERMUTE, tei0, tei0
              vmovapd        ter0, (%rax)
              vmovupd        ter2, (%rax, %rdi, 8)
              vmovapd        ter3, (%rbx)
              vmovupd        ter1, (%rbx, %rdi, 8)
              vmovapd        tei1, (%rcx)
              vmovupd        tei3, (%rcx, %rdi, 8)
              vmovapd        tei2, (%rdx)
              vmovupd        tei0, (%rdx, %rdi, 8)
              xor            k2, k2
radix4_begin: add            $4, k2
              mov            N2, %rax
              shr            %rax
              cmp            k2, %rax
              jle            radix4_end
              vmovapd        (Wr, k2, 8), cos1
              vmovapd        (Wi, k2, 8), sin1
              vmulpd         sin1, sin1, cos2
              vmulpd         cos1, sin1, sin2
              vfmsub231pd    cos1, cos1, cos2
              vfmadd231pd    sin1, cos1, sin2
              mov            N2, %rdi
              sub            k2, %rdi
              sub            k2, %rdi
              sub            $3, %rdi
              lea            (X, k2, 8), %rax
              lea            (%rax, N2, 8), %rbx
              lea            (%rbx, N2, 8), %rcx
              lea            (%rcx, N2, 8), %rdx
              vmovapd        (%rax), er0
              vmovapd        (%rbx), er1
              vmovapd        (%rcx), er2
              vmovapd        (%rdx), er3
              vpermpd        $27, (%rax, %rdi, 8), ei0
              vpermpd        $27, (%rbx, %rdi, 8), ei1
              vpermpd        $27, (%rcx, %rdi, 8), ei2
              vpermpd        $27, (%rdx, %rdi, 8), ei3
              call           butter_4
              vpermpd        $27, er2, er2
              vpermpd        $27, er1, er1
              vpermpd        $27, ei3, ei3
              vpermpd        $27, ei0, ei0
              vmovapd        er0, (%rax)
              vmovupd        er2, (%rax, %rdi, 8)
              vmovapd        er3, (%rbx)
              vmovupd        er1, (%rbx, %rdi, 8)
              vmovapd        ei1, (%rcx)
              vmovupd        ei3, (%rcx, %rdi, 8)
              vmovapd        ei2, (%rdx)
              vmovupd        ei0, (%rdx, %rdi, 8)
              jmp            radix4_begin
radix4_end:   ret

butter_4:     vmovapd        er0, tr
              vmovapd        ei0, ti
              vfmadd231pd    sin2, ei2, tr
              vfnmadd231pd   sin2, er2, ti
              vfnmadd132pd   cos2, tr, er2
              vfnmadd132pd   cos2, ti, ei2
              vfmsub132pd    TWO, er2, er0
              vfmsub132pd    TWO, ei2, ei0
              vmovapd        er1, tr
              vmovapd        ei1, ti
              vfmadd231pd    sin2, ei3, tr
              vfnmadd231pd   sin2, er3, ti
              vfnmadd132pd   cos2, tr, er3
              vfnmadd132pd   cos2, ti, ei3
              vfmsub132pd    TWO, er3, er1
              vfmsub132pd    TWO, ei3, ei1
              vmovapd        er0, tr
              vmovapd        ei0, ti
              vfmadd231pd    sin1, ei1, tr
              vfnmadd231pd   sin1, er1, ti
              vfnmadd132pd   cos1, tr, er1
              vfmsub132pd    cos1, ti, ei1
              vfmsub132pd    TWO, er1, er0
              vfmadd132pd    TWO, ei1, ei0
              vmovapd        er2, tr
              vmovapd        ei2, ti
              vfmadd231pd    cos1, ei3, tr
              vfnmadd231pd   cos1, er3, ti
              vfmadd132pd    sin1, tr, er3
              vfmadd132pd    sin1, ti, ei3
              vfmsub132pd    TWO, er3, er2
              vfnmadd132pd   TWO, ei3, ei2
              ret

butter_2:     vmovapd        ter0, ttr
              vmovapd        tei0, tti
              vfmadd231pd    tsin2, tei2, ttr
              vfnmadd231pd   tsin2, ter2, tti
              vfnmadd132pd   tcos2, ttr, ter2
              vfnmadd132pd   tcos2, tti, tei2
              vfmsub132pd    TWO, ter2, ter0
              vfmsub132pd    TWO, tei2, tei0
              vmovapd        ter1, ttr
              vmovapd        tei1, tti
              vfmadd231pd    tsin2, tei3, ttr
              vfnmadd231pd   tsin2, ter3, tti
              vfnmadd132pd   tcos2, ttr, ter3
              vfnmadd132pd   tcos2, tti, tei3
              vfmsub132pd    TWO, ter3, ter1
              vfmsub132pd    TWO, tei3, tei1
              vmovapd        ter0, ttr
              vmovapd        tei0, tti
              vfmadd231pd    tsin1, tei1, ttr
              vfnmadd231pd   tsin1, ter1, tti
              vfnmadd132pd   tcos1, ttr, ter1
              vfmsub132pd    tcos1, tti, tei1
              vfmsub132pd    TWO, ter1, ter0
              vfmadd132pd    TWO, tei1, tei0
              vmovapd        ter2, ttr
              vmovapd        tei2, tti
              vfmadd231pd    tcos1, tei3, ttr
              vfnmadd231pd   tcos1, ter3, tti
              vfmadd132pd    tsin1, ttr, ter3
              vfmadd132pd    tsin1, tti, tei3
              vfmsub132pd    TWO, ter3, ter2
              vfnmadd132pd   TWO, tei3, tei2
              ret

butter_1:     vmovq          ter0, ttr
              vmovq          tei0, tti
              vfmadd231sd    tsin2, tei2, ttr
              vfnmadd231sd   tsin2, ter2, tti
              vfnmadd132sd   tcos2, ttr, ter2
              vfnmadd132sd   tcos2, tti, tei2
              vfmsub132sd    TWO, ter2, ter0
              vfmsub132sd    TWO, tei2, tei0
              vmovq          ter1, ttr
              vmovq          tei1, tti
              vfmadd231sd    tsin2, tei3, ttr
              vfnmadd231sd   tsin2, ter3, tti
              vfnmadd132sd   tcos2, ttr, ter3
              vfnmadd132sd   tcos2, tti, tei3
              vfmsub132sd    TWO, ter3, ter1
              vfmsub132sd    TWO, tei3, tei1
              vmovq          ter0, ttr
              vmovq          tei0, tti
              vfmadd231sd    tsin1, tei1, ttr
              vfnmadd231sd   tsin1, ter1, tti
              vfnmadd132sd   tcos1, ttr, ter1
              vfmsub132sd    tcos1, tti, tei1
              vfmsub132sd    TWO, ter1, ter0
              vfmadd132sd    TWO, tei1, tei0
              vmovq          ter2, ttr
              vmovq          tei2, tti
              vfmadd231sd    tcos1, tei3, ttr
              vfnmadd231sd   tcos1, ter3, tti
              vfmadd132sd    tsin1, ttr, ter3
              vfmadd132sd    tsin1, tti, tei3
              vfmsub132sd    TWO, ter3, ter2
              vfnmadd132sd   TWO, tei3, tei2
              ret

butter0_4:    vaddpd         er2, er0, ei0
              vsubpd         er2, er0, ei2
              vaddpd         er3, er1, ei1
              vsubpd         er1, er3, er3
              vaddpd         ei1, ei0, er0
              vsubpd         ei1, ei0, er2
              vmovapd        ei2, er1
              ret

butter0_1:    vaddsd         ter2, ter0, tei0
              vsubsd         ter2, ter0, tei2
              vaddsd         ter3, ter1, tei1
              vsubsd         ter1, ter3, ter3
              vaddsd         tei1, tei0, ter0
              vsubsd         tei1, tei0, ter2
              vmovq          tei2, ter1
              ret

butterNy_1:   vaddsd         ter3, ter1, tei0
              vsubsd         ter3, ter1, tei2
              vmulsd         TW45, tei2, tei1
              vmulsd         TW45, tei0, tei3
              vmovq          ter0, tei0
              vmovq          ter2, tei2
              vaddsd         tei1, tei0, ter0
              vaddsd         tei3, tei2, ter3
              vsubsd         tei1, tei0, ter1
              vsubsd         tei3, tei2, ter2
              vmulpd         NONE, ter3, ter3
              ret


              .align         32
TWO:          .double        2.0
              .double        2.0
              .double        2.0
              .double        2.0
NONE:         .double        -1.0
              .double        -1.0
              .double        -1.0
              .double        -1.0
TW45:         .double        .7071067812
              .double        .7071067812
              .double        .7071067812
              .double        .7071067812
PERMUTE:      .byte          8
              .byte          9
              .byte          10
              .byte          11
              .byte          12
              .byte          13
              .byte          14
              .byte          15
              .byte          0
              .byte          1
              .byte          2
              .byte          3
              .byte          4
              .byte          5
              .byte          6
              .byte          7

