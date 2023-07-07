#define       X              %r15
#define       C              %r14
#define       S              %r13
#define       NLO            %r12
#define       N              %r11
#define       N2             %r10
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


              .global        dit_nr_recur

              .text

dit_nr_recur: push           %r12
              push           %r13
              push           %r14
              push           %r15
              push           %rbx
              mov            %rdi, X
              mov            %rsi, N
              mov            N, %rdi
              push           %r11
              call           get_twiddles
              mov            %rax, C
              mov            %rdx, S
              pop            %r11
              mov            $1, N2
              xor            k2rev, k2rev
              vmovapd        TW45, tw45
              vmovapd        TWO, two
              call           recurse
              pop            %rbx
              pop            %r15
              pop            %r14
              pop            %r13
              pop            %r12
              ret
recurse:      mov            N, NLO
              shr            logN1, NLO
              test           $1, k2rev
              jnz            rcalls
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
              add            $4, ilo
              cmp            ilo, NLO
              jne            keq0_loop
              cmp            $1, N2
              jle            finish
              xor            ilo, ilo
kny_loop:     vmovapd        (%rax, N, 8), ei0
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
              vsubpd         er1, er0, ei2
              vsubpd         er3, er2, ei1
              vmulpd         NONE, ei3, ei3
              vmovapd        ei0, (%rax, N, 8)
              vmovapd        ei2, (%rbx, N, 8)
              vmovapd        ei1, (%rcx, N, 8)
              vmovapd        ei3, (%rdx, N, 8)
              jmp            finish
krest:        vbroadcastsd   (C, k2rev, 8), cos1
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
              vmovapd        (%rax, N, 8), ei0
              vmovapd        (%rax, N, 8), ei1
              vmovapd        (%rax, N, 8), ei2
              vmovapd        (%rax, N, 8), ei3
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
              vmovapd        ei0, (%rbx)
              vmovapd        er3, (%rcx)
              vmovapd        ei3, (%rdx)
              vmovapd        er2, (%rax, N, 8)
              vmovapd        ei2, (%rbx, N, 8)
              vmovapd        er1, (%rcx, N, 8)
              vmovapd        ei1, (%rdx, N, 8)
              add            $4, ilo
              cmp            ilo, NLO
              jne            krest_loop
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
              mov            $3, %rcx
              mov            $1, %rax
              cmp            $1, N2
              cmove          %rax, %rcx
scalarlp:     imul           $8, %rcx, %rax
              vmovq          (C, %rcx, 8), tcos1
              vmovq          (S, %rcx, 8), tsin1
              vmulsd         tsin1, tsin1, tcos2
              vmulsd         tcos1, tsin1, tsin2
              vfmsub231sd    tcos1, tcos1, tcos2
              vfmadd231sd    tsin1, tcos1, tsin2
              vmovq          0(X, %rax, 8), ur0
              vmovq          8(X, %rax, 8), ur1
              vmovq          16(X, %rax, 8), ur2
              vmovq          24(X, %rax, 8), ur3
              vmovq          32(X, %rax, 8), ui0
              vmovq          40(X, %rax, 8), ui1
              vmovq          48(X, %rax, 8), ui2
              vmovq          56(X, %rax, 8), ui3
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
              vmovq          0(X, %rax, 8), ur0
              vmovq          8(X, %rax, 8), ui0
              vmovq          16(X, %rax, 8), ur3
              vmovq          24(X, %rax, 8), ui3
              vmovq          32(X, %rax, 8), ur2
              vmovq          40(X, %rax, 8), ui2
              vmovq          48(X, %rax, 8), ur1
              vmovq          56(X, %rax, 8), ui1
              dec            %rcx
              jnz            scalarlp
              jmp            done
finalk2:      shl            logN1, k2rev
              vmovupd        (S, k2rev, 8), sin1
              vmovupd        (C, k2rev, 8), cos1
              vmulpd         sin1, sin1, cos2
              vmulpd         cos1, sin1, sin2
              vfmsub231pd    cos1, cos1, cos2
              vfmadd231pd    sin1, cos1, sin2
              lea            (X, N, 8), %rax
              vmovapd        none, %ymm8
              vmovapd        TINDICES, %xmm9
              vmovapd        %ymm8, %ymm10
              vmovapd        %ymm8, %ymm11
              vmovapd        %ymm8, %ymm12
              vmovapd        %ymm8, %ymm13
              vgatherdpd     %ymm10, (X, %xmm9, 8), er0
              vgatherdpd     %ymm11, 8(X, %xmm9, 8), er1
              vgatherdpd     %ymm12, 16(X, %xmm9, 8), er2
              vgatherdpd     %ymm13, 24(X, %xmm9, 8), er3
              vmovapd        %ymm8, %ymm10
              vmovapd        %ymm8, %ymm11
              vmovapd        %ymm8, %ymm12
              vmovapd        %ymm8, %ymm13
              vgatherdpd     %ymm10, (%rax, %xmm9, 8), ei0
              vgatherdpd     %ymm11, 8(%rax, %xmm9, 8), ei1
              vgatherdpd     %ymm12, 16(%rax, %xmm9, 8), ei2
              vgatherdpd     %ymm13, 24(%rax, %xmm9, 8), ei3
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
              vshufpd        $0, er3, er0, %ymm4
              vshufpd        $15, er3, er0, %ymm5
              vshufpd        $0, er1, er2, %ymm6
              vshufpd        $15, er1, er2, %ymm7
              vinsertf128    $1, %xmm6, %ymm4, er0
              vinsertf128    $1, %xmm7, %ymm5, er2
              vpermpd        $78, %ymm4, %ymm4
              vpermpd        $78, %ymm5, %ymm5
              vinsertf128    $0, %xmm4, %ymm6, er1
              vinsertf128    $0, %xmm5, %ymm7, er3
              vshufpd        $0, ei3, ei0, %ymm4
              vshufpd        $15, ei3, ei0, %ymm5
              vshufpd        $0, ei1, ei2, %ymm6
              vshufpd        $15, ei1, ei2, %ymm7
              vinsertf128    $1, %xmm6, %ymm4, ei0
              vinsertf128    $1, %xmm7, %ymm5, ei2
              vpermpd        $78, %ymm4, %ymm4
              vpermpd        $78, %ymm5, %ymm5
              vinsertf128    $0, %xmm4, %ymm6, ei1
              vinsertf128    $0, %xmm5, %ymm7, ei3
              vmovapd        er0, (X)
              vmovapd        er2, 32(X)
              vmovapd        er1, 64(X)
              vmovapd        er3, 96(X)
              vmovapd        ei0, (%rax)
              vmovapd        ei2, 32(%rax)
              vmovapd        ei1, 64(%rax)
              vmovapd        ei3, 96(%rax)
              jmp            done
rcalls:       cmp            $4, NLO
              jle            done
              push           k2rev
              push           X
              shl            logN1, k2rev
              push           k2rev
              push           X
              shr            logN1, N
              inc            k2rev
              lea            (X, N, 8), X
              push           k2rev
              push           X
              inc            k2rev
              lea            (X, N, 8), X
              push           k2rev
              push           X
              inc            k2rev
              lea            (X, N, 8), X
              shl            logN1, N2
              call           recurse
              shr            logN1, N2
              pop            X
              pop            k2rev
              shl            logN1, N2
              call           recurse
              shr            logN1, N2
              pop            X
              pop            k2rev
              shl            logN1, N2
              call           recurse
              shr            logN1, N2
              pop            X
              pop            k2rev
              shl            logN1, N2
              call           recurse
              pop            X
              pop            k2rev
              shr            logN1, N2
              shl            logN1, N
done:         ret

              .align         32
TW45:         .double        0.70710678118654752440
              .double        0.70710678118654752440
              .double        0.70710678118654752440
              .double        0.70710678118654752440
NONE:         .double        -1.0
              .double        -1.0
              .double        -1.0
              .double        -1.0
TWO:          .double        2.0
              .double        2.0
              .double        2.0
              .double        2.0

              .align         16
TINDICES:     .long          0
              .long          4
              .long          8
              .long          12

