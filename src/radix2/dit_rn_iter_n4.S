#include      "ymm.h"
#define       X              %r15
#define       Y              %r14
#define       N1N2           %r13
#define       N2             %r12
#define       nhi            %r11
#define       k2             %r10
#define       Wr             %r9
#define       Wi             %r8
#define       STACK_SIZE     $8
#define       NHI            -8(%rbp)

              .global        dit_rn_iter_n4


              .text
dit_rn_iter_n4:
              push           %rbp
              mov            %rsp, %rbp
              sub            STACK_SIZE, %rsp
              push           %rbx
              push           %r12
              push           %r13
              push           %r14
              push           %r15
              mov            %rdi, X
              mov            %rsi, Y
              mov            %rdx, NHI
              mov            %rcx, N2
              mov            N2, %rdi
              shl            $2, %rdi
              mov            %rdi, N1N2
              call           twiddles_fwd_complex
              mov            %rax, Wr
              mov            %rdx, Wi
              xor            nhi, nhi
              imul           $2, N2, %rcx
              vmovapd        TWO, two
              imul           $3, N2, %rdx
nhi_loop:     xor            k2, k2
k2_loop:      vmovapd        (Wr, k2, 8), cos1
              vmovapd        (Wi, k2, 8), sin1
              vmulpd         sin1, sin1, cos2
              vmulpd         cos1, sin1, sin2
              vfmsub231pd    cos1, cos1, cos2
              mov            nhi, %rax
              imul           N1N2, %rax
              lea            (X, %rax, 8), %rsi
              lea            (Y, %rax, 8), %rdi
              lea            (%rsi, k2, 8), %rsi
              lea            (%rdi, k2, 8), %rdi
              vmovapd        (%rsi), er0
              vmovapd        (%rsi, N1N2, 8), er2
              vmovapd        (%rsi, %rcx, 8), er1
              vmovapd        (%rsi, %rdx, 8), er3
              vmovapd        (%rdi), ei0
              vmovapd        (%rdi, N1N2, 8), ei2
              vmovapd        (%rdi, %rcx, 8), ei1
              vmovapd        (%rdi, %rdx, 8), ei3
              call           tw_butterfly4
              vmovapd        er0, (%rsi)
              vmovapd        er3, (%rsi, N1N2, 8)
              vmovapd        er1, (%rsi, %rcx, 8)
              vmovapd        er2, (%rsi, %rdx, 8)
              vmovapd        ei0, (%rdi)
              vmovapd        ei3, (%rdi, N1N2, 8)
              vmovapd        ei1, (%rdi, %rcx, 8)
              vmovapd        ei2, (%rdi, %rdx, 8)
              add            $4, k2
              cmp            k2, N2
              jne            k2_loop
              inc            nhi
              cmp            nhi, NHI
              jne            nhi_loop
              pop            %r15
              pop            %r14
              pop            %r13
              pop            %r12
              pop            %rbx
              mov            %rbp, %rsp
              pop            %rbp
              ret
              .align         32
TWO:          .double        2.0
              .double        2.0
              .double        2.0
              .double        2.0
