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
#define       none           %ymm15
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
#define       tnone          %xmm15

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
              call           next_level
              mov            X, %rdi
              mov            N, %rsi
              call           scramble
              pop            %rbx
              pop            %r15
              pop            %r14
              pop            %r13
              pop            %r12
              ret
next_level:   bsf            N, %rax
              bt             $0, %rax
              jc             radix8
              mov            N, NLO
              shr            $2, NLO
              cmp            $1, k2rev
              je             next_levels
              test           k2rev, k2rev
              jz             k2eq0
              jmp            k2_interior
k2eq0:        xor            ilo, ilo
k2eq0_loop:   lea            (X, ilo, 8), %rax
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
              jne            k2eq0_loop
              bt             $0, N2
              jc             skip_k2eqNy
              vmovapd        TW45, tw45
              vmovapd        NONE, none
              xor            ilo, ilo
k2eqNy_loop:  lea            (X, ilo, 8), %rax
              lea            (%rax, N, 8), %rax
              lea            (%rax, NLO, 8), %rbx
              lea            (%rbx, NLO, 8), %rcx
              lea            (%rcx, NLO, 8), %rdx
              vmovapd        (%rax), ei0
              vmovapd        (%rbx), ei1
              vmovapd        (%rcx), ei2
              vmovapd        (%rdx), ei3
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
              vmulpd         none, ei3, ei3
              vmovapd        ei0, (%rax)
              vmovapd        ei2, (%rbx)
              vmovapd        ei1, (%rcx)
              vmovapd        ei3, (%rdx)
              add            $4, ilo
              cmp            ilo, NLO
              jne            k2eqNy_loop
skip_k2eqNy:  cmp            $4, NLO
              jg             next_levels
              jmp            final_level
k2_interior:  mov            k2rev, %rdi
              imul           NLO, %rdi
              shl            $2, %rdi
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
              jl             next_levels
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
              vmovapd        TWO, two
              xor            ilo, ilo
k2_int_loop:  lea            (X, ilo, 8), %rax
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
              jne            k2_int_loop
              pop            X
              cmp            $4, NLO
              jg             next_levels
              jmp            final_level
next_levels:  cmp            $4, NLO
              jle            done2
              push           k2rev
              push           X
              shr            $2, N
              shl            $2, k2rev
              add            $3, k2rev
              imul           $3, N, %rax
              lea            (X, %rax, 8), %rax
              push           k2rev
              push           %rax
              dec            k2rev
              imul           $2, N, %rax
              lea            (X, %rax, 8), %rax
              push           k2rev
              push           %rax
              dec            k2rev
              lea            (X, N, 8), %rax
              push           k2rev
              push           %rax
              dec            k2rev
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
              pop            X
              pop            k2rev
              shr            $2, N2
              shl            $2, N
done2:        ret
final_level:  test           k2rev, k2rev
              jnz            final_k2_simd
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
              vmulsd         TW45, ur2, ur1
              vmulsd         TW45, ur0, ur3
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
              mov            $2, %rsi
              mov            $1, %rdi
              vmovapd        TWO, ttwo
              call           final_k2_sclr
              bt             $0, N2
              jc             done1
              mov            $4, %rsi
              mov            $3, %rdi
              call           final_k2_sclr
              mov            $6, %rsi
              mov            $-1, %rdi
              call           final_k2_sclr
done1:        ret
final_k2_sclr:vmovq          (C, %rsi, 8), tcos1
              vmovq          (S, %rsi, 8), tsin1
              shl            $2, %rsi
              shl            $2, %rdi
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
final_k2_simd:mov            k2rev, %rdi
              imul           NLO, %rdi
              shl            $2, %rdi
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
              vunpcklpd      %ymm2, %ymm0, %ymm8
              vunpcklpd      %ymm3, %ymm1, %ymm12
              vunpckhpd      %ymm2, %ymm0, %ymm9
              vunpckhpd      %ymm3, %ymm1, %ymm13
              vunpcklpd      %ymm6, %ymm4, %ymm10
              vunpcklpd      %ymm7, %ymm5, %ymm14
              vunpckhpd      %ymm6, %ymm4, %ymm11
              vunpckhpd      %ymm7, %ymm5, %ymm15
              vperm2f128     $0x20, %ymm10, %ymm8, %ymm0
              vperm2f128     $0x20, %ymm14, %ymm12, %ymm4
              vperm2f128     $0x20, %ymm11, %ymm9, %ymm1
              vperm2f128     $0x20, %ymm15, %ymm13, %ymm5
              vperm2f128     $0x31, %ymm10, %ymm8, %ymm2
              vperm2f128     $0x31, %ymm14, %ymm12, %ymm6
              vperm2f128     $0x31, %ymm11, %ymm9, %ymm3
              vperm2f128     $0x31, %ymm15, %ymm13, %ymm7
              vmovapd        TWO, two
              shl            $2, k2rev
              vpermpd        $120, (S, k2rev, 8), sin1
              vpermpd        $120, (C, k2rev, 8), cos1
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
              vunpcklpd      %ymm5, %ymm0, %ymm8
              vunpcklpd      %ymm7, %ymm2, %ymm12
              vunpckhpd      %ymm5, %ymm0, %ymm9
              vunpckhpd      %ymm7, %ymm2, %ymm13
              vunpcklpd      %ymm6, %ymm3, %ymm10
              vunpcklpd      %ymm4, %ymm1, %ymm14
              vunpckhpd      %ymm6, %ymm3, %ymm11
              vunpckhpd      %ymm4, %ymm1, %ymm15
              vperm2f128     $0x20, %ymm10, %ymm8, %ymm0
              vperm2f128     $0x20, %ymm14, %ymm12, %ymm1
              vperm2f128     $0x20, %ymm11, %ymm9, %ymm2
              vperm2f128     $0x20, %ymm15, %ymm13, %ymm3
              vperm2f128     $0x31, %ymm10, %ymm8, %ymm4
              vperm2f128     $0x31, %ymm14, %ymm12, %ymm5
              vperm2f128     $0x31, %ymm11, %ymm9, %ymm6
              vperm2f128     $0x31, %ymm15, %ymm13, %ymm7
              vmovapd        %ymm0, (X)
              vmovapd        %ymm1, 32(X)
              vmovapd        %ymm2, 64(X)
              vmovapd        %ymm3, 96(X)
              vmovapd        %ymm4, (X, %rdi, 8)
              vmovapd        %ymm5, 32(X, %rdi, 8)
              vmovapd        %ymm6, 64(X, %rdi, 8)
              vmovapd        %ymm7, 96(X, %rdi, 8)
              ret
radix8:       mov            N, NLO
              shr            $3, NLO
              xor            ilo, ilo
              vmovapd        NONE, %ymm15
              vmovapd        TW45, %ymm14
              mov            NLO, %rdi
              shl            $3, %rdi
radix8_loop:  lea            (X, ilo, 8), %rax
              lea            (%rax, %rdi, 2), %rbx
              lea            (%rbx, %rdi, 2), %rcx
              lea            (%rcx, %rdi, 2), %rdx
              vmovapd        (%rax), %ymm0
              vmovapd        (%rax, %rdi), %ymm1
              vmovapd        (%rbx), %ymm2
              vmovapd        (%rbx, %rdi), %ymm3
              vmovapd        (%rcx), %ymm4
              vmovapd        (%rcx, %rdi), %ymm5
              vmovapd        (%rdx), %ymm6
              vmovapd        (%rdx, %rdi, 1), %ymm7
	          vaddpd         %ymm0, %ymm4, %ymm8
	          vsubpd         %ymm4, %ymm0, %ymm0
	          vaddpd         %ymm2, %ymm6, %ymm4
	          vsubpd         %ymm6, %ymm2, %ymm2
	          vaddpd         %ymm8, %ymm4, %ymm6
	          vmovapd        %ymm2, %ymm9
	          vsubpd         %ymm4, %ymm8, %ymm2
	          vaddpd         %ymm1, %ymm5, %ymm8
	          vsubpd         %ymm5, %ymm1, %ymm1
	          vaddpd         %ymm7, %ymm3, %ymm4
	          vsubpd         %ymm3, %ymm7, %ymm7
	          vsubpd         %ymm4, %ymm8, %ymm3
	          vaddpd         %ymm8, %ymm4, %ymm8
	          vmovapd        %ymm0, %ymm4
	          vaddpd         %ymm6, %ymm8, %ymm0
	          vmovapd        %ymm6, %ymm5
	          vmulpd         %ymm15, %ymm3, %ymm6
	          vmovapd        %ymm4, %ymm3
	          vsubpd         %ymm8, %ymm5, %ymm4
	          vmulpd         %ymm14, %ymm7, %ymm7
	          vmulpd         %ymm14, %ymm1, %ymm1
	          vaddpd         %ymm1, %ymm7, %ymm8
	          vmovapd        %ymm1, %ymm5
	          vaddpd         %ymm3, %ymm8,%ymm1
	          vsubpd         %ymm5, %ymm7, %ymm7
	          vmovapd        %ymm7, %ymm5
	          vsubpd         %ymm9, %ymm5, %ymm7
	          vmovapd        %ymm3, %ymm10
	          vsubpd         %ymm8, %ymm10, %ymm3
	          vmovapd        %ymm5, %ymm8
	          vaddpd         %ymm8, %ymm9, %ymm5
              vmovapd        %ymm0, (%rax)
              vmovapd        %ymm4, (%rax, %rdi)
              vmovapd        %ymm2, (%rbx)
              vmovapd        %ymm6, (%rbx, %rdi)
              vmovapd        %ymm1, (%rcx)
              vmovapd        %ymm5, (%rcx, %rdi)
              vmovapd        %ymm3, (%rdx)
              vmovapd        %ymm7, (%rdx, %rdi)
              add            $4, ilo
              cmp            ilo, NLO
              jg             radix8_loop
              cmp            $4, NLO
              jle            done3
              push           k2rev
              push           X
              shr            $3, N
              shl            $3, k2rev
              add            $7, k2rev
              imul           $7, N, %rax
              lea            (X, %rax, 8), %rax
              push           k2rev
              push           %rax
              dec            k2rev
              imul           $6, N, %rax
              lea            (X, %rax, 8), %rax
              push           k2rev
              push           %rax
              dec            k2rev
              imul           $5, N, %rax
              lea            (X, %rax, 8), %rax
              push           k2rev
              push           %rax
              dec            k2rev
              imul           $4, N, %rax
              lea            (X, %rax, 8), %rax
              push           k2rev
              push           %rax
              dec            k2rev
              imul           $3, N, %rax
              lea            (X, %rax, 8), %rax
              push           k2rev
              push           %rax
              dec            k2rev
              imul           $2, N, %rax
              lea            (X, %rax, 8), %rax
              push           k2rev
              push           %rax
              dec            k2rev
              lea            (X, N, 8), %rax
              push           k2rev
              push           %rax
              dec            k2rev
              shl            $3, N2
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
              pop            X
              pop            k2rev
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
              pop            X
              pop            k2rev
              shr            $3, N2
              shl            $3, N
done3:        ret


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
