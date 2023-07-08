#define       X              %r15
#define       C              %r14
#define       S              %r13
#define       NLO            %r12
#define       N              %r11
#define       pos            %r10
#define       k2rev          %r9
#define       ilo            %r8
#define       er0            %ymm0
#define       er1            %ymm1
#define       er2            %ymm2
#define       er3            %ymm3
#define       ei0            %ymm4
#define       ei1            %ymm5
#define       ei2            %ymm6
#define       ei3            %ymm7
#define       cos1           %ymm8
#define       sin1           %ymm9
#define       cos2           %ymm10
#define       sin2           %ymm11
#define       tr             %ymm12
#define       ti             %ymm13
#define       tw45           %ymm14
#define       two            %ymm15
#define       ur0            %xmm0
#define       ur1            %xmm1
#define       ur2            %xmm2
#define       ur3            %xmm3
#define       ui0            %xmm4
#define       ui1            %xmm5
#define       ui2            %xmm6
#define       ui3            %xmm7
#define       tcos1          %xmm8
#define       tsin1          %xmm9
#define       tcos2          %xmm10
#define       tsin2          %xmm11
#define       ttr            %xmm12
#define       tti            %xmm13
#define       ttw45          %xmm14
#define       ttwo           %xmm15
#define       N1             $4
#define       logN1          $2
#define       STACK_SIZE     $16
#define       N2             -8(%rbp)


              .global        dit_nr_recur

              .text

dit_nr_recur: push           %rbp
              mov            %rsp, %rbp
              sub            STACK_SIZE, %rsp
              push           %r12
              push           %r13
              push           %r14
              push           %r15
              push           %rbx
              mov            %rdi, X
              mov            %rsi, N
              mov            $1, %rax
              mov            %rax, N2
              mov            N, %rdi
              push           %r11
              call           get_twiddles
              mov            %rax, C
              mov            %rdx, S
              pop            %r11
              xor            pos, pos
              xor            k2rev, k2rev
              vmovapd        TW45, tw45
              vmovapd        TWO, two
              call           recurse
              mov            X, %rdi
              mov            N, %rsi
              push           %r10
              push           %r11
              call           scramble
              pop            %r11
              pop            %r10
              pop            %rbx
              pop            %r15
              pop            %r14
              pop            %r13
              pop            %r12
              mov            %rbp, %rsp
              pop            %rbp
              ret
recurse:      mov            N, NLO
              shr            logN1, NLO
              cmp            $1, k2rev
              je             rcalls
              cmp            $0, k2rev
              je             keq0
              jmp            krest
keq0:         xor            ilo, ilo
keq0_loop:    lea            (X, ilo, 8), %rax
              lea            (%rax, NLO, 8), %rbx
              lea            (%rbx, NLO, 8), %rcx
              lea            (%rcx, NLO, 8), %rdx
              vmovapd        (%rax), er0
              vmovapd        (%rbx), er1
              vmovapd        (%rcx), er2
              vmovapd        (%rdx), er3
              vaddpd         er2, er0, ei0
              vsubpd         er2, er0, ei2
              vaddpd         er3, er1, ei1
              vsubpd         er1, er3, ei3
              vaddpd         ei1, ei0, er0
              vsubpd         ei1, ei0, er2
              vmovapd        er0, (%rax)
              vmovapd        er2, (%rbx)
              vmovapd        ei2, (%rcx)
              vmovapd        ei3, (%rdx)
              cmp            $1, N2
              jle            skip_kny_loop
              vmovapd        (%rax, N, 8), ei0
              vmovapd        (%rbx, N, 8), ei1
              vmovapd        (%rcx, N, 8), ei2
              vmovapd        (%rdx, N, 8), ei3
              vaddpd         ei3, ei1, er0
              vsubpd         ei3, ei1, er2
              vmulpd         tw45, er2, er1
              vmulpd         tw45, er0, er3
              vmovapd        ei0, er0
              vmovapd        ei2, er2
              vaddpd         er1, er0, ei0
              vaddpd         er3, er2, ei3
              vsubpd         er1, er0, ei1
              vsubpd         er3, er2, ei2
              vmulpd         NONE, ei3, ei3
              vmovapd        ei0, (%rax, N, 8)
              vmovapd        ei2, (%rbx, N, 8)
              vmovapd        ei1, (%rcx, N, 8)
              vmovapd        ei3, (%rdx, N, 8)
skip_kny_loop:add            $4, ilo
              cmp            ilo, NLO
              jne            keq0_loop
              jmp            finish
krest:        mov            pos, %rdi
              mov            %rdi, %rsi
              bsr            %rdi, %rcx
              mov            $3, %rdi
              shl            %rcx, %rdi
              dec            %rdi
              sub            %rsi, %rdi
              sub            %rsi, %rdi
              mov            N, %rax
              dec            %rax
              sub            %rax, %rdi
              cmp            $0, %rdi
              jl             rcalls
              push           X
              lea            (X, %rdi, 8), %rax
              mov            %rdi, %rdx
              neg            %rdx
              test           $1, k2rev
              cmovnz         %rax, X
              cmovnz         %rdx, %rdi
              vbroadcastsd   (C, k2rev, 8), cos1
              vbroadcastsd   (S, k2rev, 8), sin1
              vmulpd         sin1, sin1, cos2
              vmulpd         cos1, sin1, sin2
              vfmsub231pd    cos1, cos1, cos2
              vfmadd231pd    sin1, cos1, sin2
              xor            ilo, ilo
krest_loop:   lea            (X, ilo, 8), %rax
              lea            (%rax, NLO, 8), %rbx
              lea            (%rbx, NLO, 8), %rcx
              lea            (%rcx, NLO, 8), %rdx
              vmovapd        (%rax), er0
              vmovapd        (%rbx), er1
              vmovapd        (%rcx), er2
              vmovapd        (%rdx), er3
              vmovapd        (%rax, %rdi, 8), ei0
              vmovapd        (%rbx, %rdi, 8), ei1
              vmovapd        (%rcx, %rdi, 8), ei2
              vmovapd        (%rdx, %rdi, 8), ei3
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
              vfmsub132pd    cos1, ti, ei1
              vfmsub132pd    two, er1, er0
              vfmadd132pd    two, ei1, ei0
              vmovapd        er2, tr
              vmovapd        ei2, ti
              vfmadd231pd    cos1, ei3, tr
              vfnmadd231pd   cos1, er3, ti
              vfmadd132pd    sin1, tr, er3
              vfmadd132pd    sin1, ti, ei3
              vfmsub132pd    two, er3, er2
              vfnmadd132pd   two, ei3, ei2
              vmovapd        er0, (%rax)
              vmovapd        ei1, (%rbx)
              vmovapd        er3, (%rcx)
              vmovapd        ei2, (%rdx)
              vmovapd        er2, (%rax, %rdi, 8)
              vmovapd        ei3, (%rbx, %rdi, 8)
              vmovapd        er1, (%rcx, %rdi, 8)
              vmovapd        ei0, (%rdx, %rdi, 8)
              add            $4, ilo
              cmp            ilo, NLO
              jne            krest_loop
              pop            X
finish:       cmp            $4, NLO
              jg             rcalls
              cmp            $0, k2rev
              jne            finalk2
              vmovq          (X), ur0
              vmovq          8(X), ur1
              vmovq          16(X), ur2
              vmovq          24(X), ur3
              vaddsd         ur2, ur0, ui0
              vsubsd         ur2, ur0, ui2
              vaddsd         ur3, ur1, ui1
              vsubsd         ur1, ur3, ui3
              vaddsd         ui1, ui0, ur0
              vsubsd         ui1, ui0, ur2
              vmovq          ur0, (X)
              vmovq          ur2, 8(X)
              vmovq          ui2, 16(X)
              vmovq          ui3, 24(X)
              vmovq          32(X), ui0
              vmovq          40(X), ui1
              vmovq          48(X), ui2
              vmovq          56(X), ui3
              vaddsd         ui3, ui1, ur0
              vsubsd         ui3, ui1, ur2
              vmulsd         ttw45, ur2, ur1
              vmulsd         ttw45, ur0, ur3
              vmovq          ui0, ur0
              vmovq          ui2, ur2
              vaddsd         ur1, ur0, ui0
              vaddsd         ur3, ur2, ui3
              vsubsd         ur1, ur0, ui1
              vsubsd         ur3, ur2, ui2
              vmulsd         NONE, ui3, ui3
              vmovq          ui0, 32(X)
              vmovq          ui2, 40(X)
              vmovq          ui1, 48(X)
              vmovq          ui3, 56(X)
scalarlp:     mov            $2, %rsi
              mov            $1, %rdi
              call           scalar_int
              cmp            $1, N2
              je             done
              mov            $4, %rsi
              mov            $3, %rdi
              call           scalar_int
              mov            $6, %rsi
              mov            $-1, %rdi
              call           scalar_int
              jmp            done
finalk2:      mov            pos, %rdi
              mov            %rdi, %rsi
              bsr            %rdi, %rcx
              mov            $3, %rdi
              shl            %rcx, %rdi
              dec            %rdi
              sub            %rsi, %rdi
              sub            %rsi, %rdi
              sub            $15, %rdi
              vmovapd        (X), %ymm0
              vmovapd        32(X), %ymm7
              vmovapd        64(X), %ymm2
              vmovapd        96(X), %ymm5
              vmovapd        (X, %rdi, 8), %ymm4
              vmovapd        32(X, %rdi, 8), %ymm3
              vmovapd        64(X, %rdi, 8), %ymm6
              vmovapd        96(X, %rdi, 8), %ymm1
              vshufpd        $0, %ymm2, %ymm0, %ymm8
              vshufpd        $0, %ymm3, %ymm1, %ymm12
              vshufpd        $15, %ymm2, %ymm0, %ymm9
              vshufpd        $15, %ymm3, %ymm1, %ymm13
              vshufpd        $0, %ymm6, %ymm4, %ymm10
              vshufpd        $0, %ymm7, %ymm5, %ymm14
              vshufpd        $15, %ymm6, %ymm4, %ymm11
              vshufpd        $15, %ymm7, %ymm5, %ymm15
              vinsertf128    $1, %xmm10, %ymm8, %ymm0
              vinsertf128    $1, %xmm14, %ymm12, %ymm4
              vinsertf128    $1, %xmm11, %ymm9, %ymm1
              vinsertf128    $1, %xmm15, %ymm13, %ymm5
              vpermpd        $78, %ymm8, %ymm8
              vpermpd        $78, %ymm12, %ymm12
              vpermpd        $78, %ymm9, %ymm9
              vpermpd        $78, %ymm13, %ymm13
              vinsertf128    $0, %xmm8, %ymm10, %ymm2
              vinsertf128    $0, %xmm12, %ymm14, %ymm6
              vinsertf128    $0, %xmm9, %ymm11, %ymm3
              vinsertf128    $0, %xmm13, %ymm15, %ymm7
              shl            logN1, k2rev
              vmovapd        EVENS, %xmm15
              vmovapd        NONE, %ymm0
              vmovapd        %ymm0, %ymm1
              and            $0xfffffffffffffffe, k2rev
              lea            (S, k2rev, 8), %rax
              lea            (C, k2rev, 8), %rbx
              vgatherdpd     %ymm0, (%rax, %xmm15, 8), sin1
              vgatherdpd     %ymm1, (%rbx, %xmm15, 8), cos1
              vmovapd        TW45, tw45
              vmovapd        TWO, two
              vmulpd         sin1, sin1, cos2
              vmulpd         cos1, sin1, sin2
              vfmsub231pd    cos1, cos1, cos2
              vfmadd231pd    sin1, cos1, sin2
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
              vfmsub132pd    cos1, ti, ei1
              vfmsub132pd    two, er1, er0
              vfmadd132pd    two, ei1, ei0
              vmovapd        er2, tr
              vmovapd        ei2, ti
              vfmadd231pd    cos1, ei3, tr
              vfnmadd231pd   cos1, er3, ti
              vfmadd132pd    sin1, tr, er3
              vfmadd132pd    sin1, ti, ei3
              vfmsub132pd    two, er3, er2
              vfnmadd132pd   two, ei3, ei2
              vpermpd        $27, %ymm4, %ymm4
              vpermpd        $27, %ymm7, %ymm7
              vpermpd        $27, %ymm1, %ymm1
              vpermpd        $27, %ymm2, %ymm2
              vshufpd        $0, %ymm5, %ymm0, %ymm8
              vshufpd        $0, %ymm7, %ymm2, %ymm12
              vshufpd        $15, %ymm5, %ymm0, %ymm9
              vshufpd        $15, %ymm7, %ymm2, %ymm13
              vshufpd        $0, %ymm6, %ymm3, %ymm10
              vshufpd        $0, %ymm4, %ymm1, %ymm14
              vshufpd        $15, %ymm6, %ymm3, %ymm11
              vshufpd        $15, %ymm4, %ymm1, %ymm15
              vinsertf128    $1, %xmm10, %ymm8, %ymm0
              vinsertf128    $1, %xmm14, %ymm12, %ymm1
              vinsertf128    $1, %xmm11, %ymm9, %ymm2
              vinsertf128    $1, %xmm15, %ymm13, %ymm3
              vpermpd        $78, %ymm8, %ymm8
              vpermpd        $78, %ymm12, %ymm12
              vpermpd        $78, %ymm9, %ymm9
              vpermpd        $78, %ymm13, %ymm13
              vinsertf128    $0, %xmm8, %ymm10, %ymm4
              vinsertf128    $0, %xmm12, %ymm14, %ymm5
              vinsertf128    $0, %xmm9, %ymm11, %ymm6
              vinsertf128    $0, %xmm13, %ymm15, %ymm7
              vmovapd        TW45, tw45
              vmovapd        TWO, two
              vmovapd        %ymm0, (X)
              vmovapd        %ymm1, 32(X)
              vmovapd        %ymm2, 64(X)
              vmovapd        %ymm3, 96(X)
              vmovapd        %ymm4, (X, %rdi, 8)
              vmovapd        %ymm5, 32(X, %rdi, 8)
              vmovapd        %ymm6, 64(X, %rdi, 8)
              vmovapd        %ymm7, 96(X, %rdi, 8)
              jmp            done
rcalls:       cmp            $4, NLO
              jle            done
              push           k2rev
              push           X
              push           pos
              shr            logN1, N
              shl            logN1, k2rev
              add            $3, k2rev
              imul           $3, N, %rax
              add            %rax, pos
              lea            (X, %rax, 8), %rax
              push           k2rev
              push           %rax
              push           pos
              dec            k2rev
              sub            N, pos
              imul           $2, N, %rax
              lea            (X, %rax, 8), %rax
              push           k2rev
              push           %rax
              push           pos
              dec            k2rev
              sub            N, pos
              lea            (X, N, 8), %rax
              push           k2rev
              push           %rax
              push           pos
              dec            k2rev
              sub            N, pos
              mov            N2, %rax
              shl            logN1, %rax
              mov            %rax, N2
              call           recurse
              pop            pos
              pop            X
              pop            k2rev
              call           recurse
              pop            pos
              pop            X
              pop            k2rev
              call           recurse
              pop            pos
              pop            X
              pop            k2rev
              call           recurse
              pop            pos
              pop            X
              pop            k2rev
              mov            N2, %rax
              shr            logN1, %rax
              mov            %rax, N2
              shl            logN1, N
done:         ret
scalar_int:   vmovq          (C, %rsi, 8), tcos1
              vmovq          (S, %rsi, 8), tsin1
              shl            logN1, %rsi
              shl            logN1, %rdi
              vmulsd         tsin1, tsin1, tcos2
              vmulsd         tcos1, tsin1, tsin2
              vfmsub231sd    tcos1, tcos1, tcos2
              vfmadd231sd    tsin1, tcos1, tsin2
              lea            (X, %rsi, 8), %rsi
              vmovq          0(%rsi), ur0
              vmovq          8(%rsi), ur1
              vmovq          16(%rsi), ur2
              vmovq          24(%rsi), ur3
              vmovq          0(%rsi, %rdi, 8), ui0
              vmovq          8(%rsi, %rdi, 8), ui1
              vmovq          16(%rsi, %rdi, 8), ui2
              vmovq          24(%rsi, %rdi, 8), ui3
              vmovq          ur0, ttr
              vmovq          ui0, tti
              vfmadd231sd    tsin2, ui2, ttr
              vfnmadd231sd   tsin2, ur2, tti
              vfnmadd132sd   tcos2, ttr, ur2
              vfnmadd132sd   tcos2, tti, ui2
              vfmsub132sd    ttwo, ur2, ur0
              vfmsub132sd    ttwo, ui2, ui0
              vmovq          ur1, ttr
              vmovq          ui1, tti
              vfmadd231sd    tsin2, ui3, ttr
              vfnmadd231sd   tsin2, ur3, tti
              vfnmadd132sd   tcos2, ttr, ur3
              vfnmadd132sd   tcos2, tti, ui3
              vfmsub132sd    ttwo, ur3, ur1
              vfmsub132sd    ttwo, ui3, ui1
              vmovq          ur0, ttr
              vmovq          ui0, tti
              vfmadd231sd    tsin1, ui1, ttr
              vfnmadd231sd   tsin1, ur1, tti
              vfnmadd132sd   tcos1, ttr, ur1
              vfmsub132sd    tcos1, tti, ui1
              vfmsub132sd    ttwo, ur1, ur0
              vfmadd132sd    ttwo, ui1, ui0
              vmovq          ur2, ttr
              vmovq          ui2, tti
              vfmadd231sd    tcos1, ui3, ttr
              vfnmadd231sd   tcos1, ur3, tti
              vfmadd132sd    tsin1, ttr, ur3
              vfmadd132sd    tsin1, tti, ui3
              vfmsub132sd    ttwo, ur3, ur2
              vfnmadd132sd   ttwo, ui3, ui2
              vmovq          ur0, 0(%rsi)
              vmovq          ui1, 8(%rsi)
              vmovq          ur3, 16(%rsi)
              vmovq          ui2, 24(%rsi)
              vmovq          ur2, 0(%rsi, %rdi, 8)
              vmovq          ui3, 8(%rsi, %rdi, 8)
              vmovq          ur1, 16(%rsi, %rdi, 8)
              vmovq          ui0, 24(%rsi, %rdi, 8)
              ret


              .align         32
TW45:         .double        0.70710678118654752440
              .double        0.70710678118654752440
              .double        0.70710678118654752440
              .double        0.70710678118654752440
NONE:         .double        -1.0
              .double        -1.0
              .double        -1.0
              .double        -1.0
ONE:          .double        1.0
              .double        1.0
              .double        1.0
              .double        1.0
ZERO:         .double        0.0
              .double        0.0
              .double        0.0
              .double        0.0
TWO:          .double        2.0
              .double        2.0
              .double        2.0
              .double        2.0

              .align         16
EVENS:        .long          0
              .long          2
              .long          4
              .long          6

