#define       k2             %r8
#define       N2o2           %r9
#define       X              %r10
#define       C              %r11
#define       S              %r12
#define       N              %r13
#define       N2             %r14
#define       NTW            %r15
#define       N1             $4
#define       logN1          $2
#define       STACK_SIZE     $16
#define       Cptrs          -8(%rbp)
#define       Sptrs          -16(%rbp)
#define       er0            %ymm0
#define       er1            %ymm1
#define       er2            %ymm2
#define       er3            %ymm3
#define       ei0            %ymm4
#define       ei1            %ymm5
#define       ei2            %ymm6
#define       ei3            %ymm7
#define       tr             %ymm8
#define       ti             %ymm9
#define       cos1           %ymm10
#define       sin1           %ymm11
#define       cos2           %ymm12
#define       sin2           %ymm13
#define       two            %ymm14
#define       none           %ymm15


              .global        dit_rn


              .data

              .align         32
NONE:         .double        -1.0
              .double        -1.0
              .double        -1.0
              .double        -1.0
TWO:          .double        2.0
              .double        2.0
              .double        2.0
              .double        2.0
SQRT1_2:      .double        0.707106781186547524400844
              .double        0.707106781186547524400844
              .double        0.707106781186547524400844
              .double        0.707106781186547524400844

              .align         16
TRANSPOSE:    .long          0
              .long          8
              .long          4
              .long          12


              .text

dit_rn:       push           %rbp
              mov            %rsp, %rbp
              sub            STACK_SIZE, %rsp
              push           %r15
              push           %r14
              push           %r13
              push           %r12
              push           %rbx
              mov            %rdi, X
              mov            %rsi, Cptrs
              mov            %rdx, Sptrs
              mov            %rcx, N
              mov            N, N2
              mov            $1, NTW
              shr            logN1, N2
              vmovapd        NONE, none
              vmovapd        TWO, two
              call           recurse
              pop            %rbx
              pop            %r12
              pop            %r13
              pop            %r14
              pop            %r15
              mov            %rbp, %rsp
              pop            %rbp
              ret
recurse:      cmp            $4, N2
              jle            N2eq4
              push           X
              mov            N2, %rdx
              neg            %rdx
              lea            (X, N, 8), %rax
              lea            (%rax, %rdx, 8), %rax
              push           %rax
              lea            (%rax, %rdx, 8), %rax
              push           %rax
              lea            (%rax, %rdx, 8), %rax
              push           %rax
              lea            (%rax, %rdx, 8), X
              shl            logN1, N2
              shl            logN1, N
              shr            logN1, NTW
              call           recurse
              pop            X
              call           recurse
              pop            X
              call           recurse
              pop            X
              call           recurse
              shl            logN1, NTW
              shr            logN1, N2
              shr            logN1, N
              pop            X
              cmp            $64, N
              jne            Ngt64
              mov            $4, %rdi
              xor            k2, k2
Neq64:        vmovapd        TRANSPOSE, %xmm8
              vmovapd        none, %ymm10
              vmovapd        none, %ymm11
              vmovapd        none, %ymm12
              vmovapd        none, %ymm13
              lea            (X, k2, 8), %rax
              vgatherdpd     %ymm10, (%rax, %xmm8, 8), er0
              vgatherdpd     %ymm11, 8(%rax, %xmm8, 8), er2
              vgatherdpd     %ymm12, 16(%rax, %xmm8, 8), er1
              vgatherdpd     %ymm13, 24(%rax, %xmm8, 8), er3
              vaddpd         er2, er0, ei0
              vsubpd         er2, er0, ei2
              vaddpd         er3, er1, ei1
              vsubpd         er1, er3, ei3
              vaddpd         ei1, ei0, er0
              vsubpd         ei1, ei0, er2
              vshufpd        $0, er2, er0, %ymm8
              vshufpd        $15, er2, er0, %ymm9
              vshufpd        $0, ei3, ei2, %ymm10
              vshufpd        $15, ei3, ei2, %ymm11
              vinsertf128    $1, %xmm10, %ymm8, er0
              vinsertf128    $1, %xmm11, %ymm9, er2
              vpermpd        $78, %ymm8, %ymm4
              vpermpd        $78, %ymm9, %ymm5
              vinsertf128    $0, %xmm8, %ymm10, er1
              vinsertf128    $0, %xmm9, %ymm11, er3
              lea            (%rax, %rdi, 8), %rbx
              lea            (%rbx, %rdi, 8), %rcx
              lea            (%rcx, %rdi, 8), %rdx
              vmovapd        er0, (%rax)
              vmovapd        er2, (%rbx)
              vmovapd        er1, (%rcx)
              vmovapd        er3, (%rdx)
              add            $4, k2
              cmp            k2, N
              jne            Neq64
              lea            (X, %rdi, 8), %rbx
              lea            (%rbx, %rdi, 8), %rcx
              lea            (%rcx, %rdi, 8), %rdx
              jmp            done
Ngt64:

done:         ret

              cmp            $4, N2
              jg             N2gt4
N2eq4:        vmovapd        TRANSPOSE, %xmm8
              vmovapd        none, %ymm10
              vmovapd        none, %ymm11
              vmovapd        none, %ymm12
              vmovapd        none, %ymm13
              vgatherdpd     %ymm10, (X, %xmm8, 8), er0
              vgatherdpd     %ymm11, 8(X, %xmm8, 8), er2
              vgatherdpd     %ymm12, 16(X, %xmm8, 8), er1
              vgatherdpd     %ymm13, 24(X, %xmm8, 8), er3
              vaddpd         er2, er0, ei0
              vsubpd         er2, er0, ei2
              vaddpd         er3, er1, ei1
              vsubpd         er1, er3, ei3
              vaddpd         ei1, ei0, er0
              vsubpd         ei1, ei0, er2
              vshufpd        $0, er2, er0, %ymm8
              vshufpd        $15, er2, er0, %ymm9
              vshufpd        $0, ei3, ei2, %ymm10
              vshufpd        $15, ei3, ei2, %ymm11
              vinsertf128    $1, %xmm10, %ymm8, er0
              vinsertf128    $1, %xmm11, %ymm9, er1
              vpermpd        $78, %ymm8, %ymm4
              vpermpd        $78, %ymm9, %ymm5
              vinsertf128    $0, %xmm8, %ymm10, er2
              vinsertf128    $0, %xmm9, %ymm11, er3
              lea            (X, N2, 8), %rbx
              lea            (%rbx, N2, 8), %rcx
              lea            (%rcx, N2, 8), %rdx
              jmp            skipN2gt4ld
N2gt4:        lea            (X, N2, 8), %rbx
              lea            (%rbx, N2, 8), %rcx
              lea            (%rcx, N2, 8), %rdx
              vmovapd        (X), er0
              vmovapd        (%rbx), er2
              vmovapd        (%rcx), er1
              vmovapd        (%rdx), er3
skipN2gt4ld:  vaddpd         er2, er0, ei0
              vsubpd         er2, er0, ei2
              vaddpd         er3, er1, ei1
              vsubpd         er1, er3, ei3
              vaddpd         ei1, ei0, er0
              vsubpd         ei1, ei0, er2
              vmovapd        er0, (X)
              vmovapd        ei2, (%rbx)
              vmovapd        er2, (%rcx)
              vmovapd        ei3, (%rdx)
              lea            (X, N2, 4), %rax
              lea            (%rax, N2, 8), %rbx
              lea            (%rbx, N2, 8), %rcx
              lea            (%rcx, N2, 8), %rdx
              vmovapd        (X), er0
              vmovapd        (%rbx), er2
              vmovapd        (%rcx), er1
              vmovapd        (%rdx), er3
              vaddpd         er3, er1, ei0
              vsubpd         er3, er1, ei2
              vmulpd         SQRT1_2, ei2, ei1
              vmulpd         SQRT1_2, ei0, ei3
              vmovapd        er0, ei0
              vmovapd        er2, ei2
              vaddpd         ei1, ei0, er0
              vaddpd         ei3, ei2, er3
              vsubpd         ei1, ei0, er2
              vsubpd         ei3, ei2, er1
              vmulpd         ei2, er3, er3
              vmovapd        er0, (X)
              vmovapd        er1, (%rbx)
              vmovapd        er2, (%rcx)
              vmovapd        er3, (%rdx)
              mov            $1, k2
              mov            N2, N2o2
              shr            N2o2
              bsf            N, %rax
              mov            Cptrs, %rdx
              mov            Sptrs, %rbx
              mov            (%rdx, %rax, 8), C
              mov            (%rbx, %rax, 8), S
k2loop:       vmovapd        (C, k2, 8), cos1
              vmovapd        (S, k2, 8), sin1
              vmulpd         sin1, sin1, cos2
              vmulpd         cos1, sin1, sin2
              vfmsub231pd    cos1, cos1, cos2
              vfmadd231pd    sin1, cos1, sin2

              add            $4, k2
              cmp            k2, N2o2
              jne            k2loop
done:         ret
