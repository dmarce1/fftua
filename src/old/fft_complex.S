#include      "ymm.h"
#define       X              %r15
#define       Y              %r14
#define       N2             %r13
#define       Wr             %r12
#define       Wi             %r11
#define       k2             %r10
#define       Wptr           %r9
#define       stack_size     %r8
#define       NHI            %rsi
#define       nhi            %rdi
#define       L1_size        $0x8000
#define       STACK_SIZE     $8
#define       N              -8(%rbp)


              .global        fft_radix2_complex


              .text

fft_radix2_complex:
              push           %rbx
              push           %r15
              push           %r14
              push           %r13
              push           %r12
              mov            %rdi, X
              mov            %rsi, Y
              mov            %rdx, N2
              shr            $2, N2
              bsf            N2, %rax
              add            $3, %rax
              shl            $4, %rax
              push           %rbp
              mov            %rsp, %rbp
              add            STACK_SIZE, %rax
              sub            %rax, %rsp
              mov            %rsp, Wptr
              mov            %rax, stack_size
              bsf            N2, %rdx
              add            $2, %rdx
              mov            $2, %rdi
twiddle_loop: mov            %rdi, %rcx
              push           %rdi
              mov            $1, %rdi
              shl            %rcx, %rdi
              push           %rdx
              push           %r8
              push           %r9
              push           %r10
              push           %r11
              call           twiddles_fwd_complex
              pop            %r11
              pop            %r10
              pop            %r9
              pop            %r8
              pop            %rcx
              pop            %rdi
              mov            %rdi, %rsi
              shl            %rsi
              mov            %rax, (Wptr, %rsi, 8)
              mov            %rdx, 8(Wptr, %rsi, 8)
              mov            %rcx, %rdx
              cmp            %rdi, %rdx
              je             twiddle_esc
              inc            %rdi
              jmp            twiddle_loop
twiddle_esc:  bsr            N2, %rax
              bt             $0, %rax
              jc             radix2
              call           next_level_r4
              jmp            cleanup
radix2:       shl            N2
              call           next_level_r2
cleanup:      add            stack_size, %rsp
              mov            %rsp, %rbp
              pop            %rbp
              pop            %r12
              pop            %r13
              pop            %r14
              pop            %r15
              pop            %rbx
              ret
next_level_r4:mov            N2, %rax
              shl            $6, %rax
              cmp            L1_size, %rax
              jl             in_l1_cache
              push           Y
              push           X
              imul           $3, N2, %rdx
              lea            (X, %rdx, 8), %rax
              lea            (Y, %rdx, 8), %rbx
              push           %rbx
              push           %rax
              imul           $2, N2, %rdx
              lea            (X, %rdx, 8), %rax
              lea            (Y, %rdx, 8), %rbx
              push           %rbx
              push           %rax
              lea            (X, N2, 8), %rax
              lea            (Y, N2, 8), %rbx
              push           %rbx
              push           %rax
              shr            $2, N2
              call           next_level_r4
              pop            X
              pop            Y
              call           next_level_r4
              pop            X
              pop            Y
              call           next_level_r4
              pop            X
              pop            Y
              call           next_level_r4
              pop            X
              pop            Y
              shl            $2, N2
              bsf            N2, %rax
              add            $2, %rax
              shl            %rax
              mov            (Wptr, %rax, 8), Wr
              mov            8(Wptr, %rax, 8), Wi
              imul           $2, N2, %rcx
              imul           $3, N2, %rdx
              vmovapd        TWO, two
              xor            k2, k2
k2_recur_loop:vmovapd        (Wr, k2, 8), cos1
              vmovapd        (Wi, k2, 8), sin1
              vmulpd         sin1, sin1, cos2
              vmulpd         cos1, sin1, sin2
              vfmsub231pd    cos1, cos1, cos2
              vfmadd231pd    sin1, cos1, sin2
              lea            (X, k2, 8), %rax
              lea            (Y, k2, 8), %rbx
              vmovapd        (%rax), er0
              vmovapd        (%rax, N2, 8), er2
              vmovapd        (%rax, %rcx, 8), er1
              vmovapd        (%rax, %rdx, 8), er3
              vmovapd        (%rbx), ei0
              vmovapd        (%rbx, N2, 8), ei2
              vmovapd        (%rbx, %rcx, 8), ei1
              vmovapd        (%rbx, %rdx, 8), ei3
              call           tw_butterfly4
              vmovapd        er0, (%rax)
              vmovapd        er3, (%rax, N2, 8)
              vmovapd        er1, (%rax, %rcx, 8)
              vmovapd        er2, (%rax, %rdx, 8)
              vmovapd        ei0, (%rbx)
              vmovapd        ei3, (%rbx, N2, 8)
              vmovapd        ei1, (%rbx, %rcx, 8)
              vmovapd        ei2, (%rbx, %rdx, 8)
              add            $4, k2
              cmp            k2, N2
              jne            k2_recur_loop
              jmp            done
in_l1_cache:  push           N2
              mov            N2, %rax
              shl            $2, %rax
              mov            %rax, N
              mov            $4, N2
              bsr            N, %rcx
              sub            $4, %rcx
              mov            $1, NHI
              shl            %rcx, NHI
              xor            nhi, nhi
nhi_loop16:   imul           $4, nhi, %rax
              imul           N2, %rax
              lea            (X, %rax, 8), %rbx
              lea            (Y, %rax, 8), %rcx
              vmovapd        (%rbx), er0
              vmovapd        32(%rbx), er1
              vmovapd        64(%rbx), er2
              vmovapd        96(%rbx), er3
              vmovapd        (%rcx), ei0
              vmovapd        32(%rcx), ei1
              vmovapd        64(%rcx), ei2
              vmovapd        96(%rcx), ei3
              call           transpose_lo
              call           butterfly4
              call           transpose_lo
              vmovapd        TWO, two
              vmovapd        tw16_1r, cos1
              vmovapd        tw16_1i, sin1
              vmovapd        tw16_2r, cos2
              vmovapd        tw16_2i, sin2
              call           tw_butterfly4
              vmovapd        er0, (%rbx)
              vmovapd        er3, 32(%rbx)
              vmovapd        er1, 64(%rbx)
              vmovapd        er2, 96(%rbx)
              vmovapd        ei0, (%rcx)
              vmovapd        ei3, 32(%rcx)
              vmovapd        ei1, 64(%rcx)
              vmovapd        ei2, 96(%rcx)
              inc            nhi
              cmp            nhi, NHI
              jne            nhi_loop16
              shl            $2, N2
              shr            $2, NHI
iter_loop:    vmovapd        TWO, two
              cmp            N, N2
              je             exit_iter
              bsf            N2, %rax
              add            $2, %rax
              shl            %rax
              mov            (Wptr, %rax, 8), Wr
              mov            8(Wptr, %rax, 8), Wi
              imul           $2, N2, %rcx
              imul           $3, N2, %rdx
              xor            k2, k2
k2_iter_loop: vmovapd        (Wr, k2, 8), cos1
              vmovapd        (Wi, k2, 8), sin1
              vmulpd         sin1, sin1, cos2
              vmulpd         cos1, sin1, sin2
              vfmsub231pd    cos1, cos1, cos2
              vfmadd231pd    sin1, cos1, sin2
              xor            nhi, nhi
nhi_loop:     imul           $4, nhi, %rbx
              imul           N2, %rbx
              add            k2, %rbx
              lea            (X, %rbx, 8), %rax
              lea            (Y, %rbx, 8), %rbx
              vmovapd        (%rax), er0
              vmovapd        (%rax, N2, 8), er2
              vmovapd        (%rax, %rcx, 8), er1
              vmovapd        (%rax, %rdx, 8), er3
              vmovapd        (%rbx), ei0
              vmovapd        (%rbx, N2, 8), ei2
              vmovapd        (%rbx, %rcx, 8), ei1
              vmovapd        (%rbx, %rdx, 8), ei3
              call           tw_butterfly4
              vmovapd        er0, (%rax)
              vmovapd        er3, (%rax, N2, 8)
              vmovapd        er1, (%rax, %rcx, 8)
              vmovapd        er2, (%rax, %rdx, 8)
              vmovapd        ei0, (%rbx)
              vmovapd        ei3, (%rbx, N2, 8)
              vmovapd        ei1, (%rbx, %rcx, 8)
              vmovapd        ei2, (%rbx, %rdx, 8)
              inc            nhi
              cmp            nhi, NHI
              jne            nhi_loop
              add            $4, k2
              cmp            k2, N2
              jne            k2_iter_loop
              shl            $2, N2
              shr            $2, NHI
              jmp            iter_loop
exit_iter:    pop            N2
done:         ret
next_level_r2:push           Y
              push           X
              lea            (X, N2, 8), %rax
              lea            (Y, N2, 8), %rbx
              push           %rbx
              push           %rax
              shr            $2, N2
              call           next_level_r4
              pop            X
              pop            Y
              call           next_level_r4
              pop            X
              pop            Y
              shl            $2, N2
              bsf            N2, %rax
              inc            %rax
              shl            %rax
              mov            (Wptr, %rax, 8), Wr
              mov            8(Wptr, %rax, 8), Wi
              vmovapd        TWO, two
              xor            k2, k2
k2_loop_r2:   vmovapd        (Wr, k2, 8), cos1
              vmovapd        (Wi, k2, 8), sin1
              lea            (X, k2, 8), %rax
              lea            (Y, k2, 8), %rbx
              vmovapd        (%rax), er0
              vmovapd        (%rax, N2, 8), er1
              vmovapd        (%rbx), ei0
              vmovapd        (%rbx, N2, 8), ei1
              vmovapd        er0, tr
              vmovapd        ei0, ti
              vfmadd231pd    sin1, ei1, tr
              vfnmadd231pd   sin1, er1, ti
              vfnmadd132pd   cos1, tr, er1
              vfnmadd132pd   cos1, ti, ei1
              vfmsub132pd    two, er1, er0
              vfmsub132pd    two, ei1, ei0
              vmovapd        er0, (%rax)
              vmovapd        er1, (%rax, N2, 8)
              vmovapd        ei0, (%rbx)
              vmovapd        ei1, (%rbx, N2, 8)
              add            $4, k2
              cmp            k2, N2
              jne            k2_loop_r2
              ret
              .align         32
TWO:          .double        2.0
              .double        2.0
              .double        2.0
              .double        2.0
tw16_1r:      .double        +1.000000000000000000
              .double        +0.923879532511286756
              .double        +0.707106781186547524
              .double        +0.382683432365089771
tw16_1i:      .double        -0.000000000000000000
              .double        -0.382683432365089771
              .double        -0.707106781186547524
              .double        -0.923879532511286756
tw16_2r:      .double        +1.000000000000000000
              .double        +0.707106781186547524
              .double        +0.000000000000000000
              .double        -0.707106781186547524
tw16_2i:      .double        +0.000000000000000000
              .double        -0.707106781186547524
              .double        -1.000000000000000000
              .double        -0.707106781186547524
