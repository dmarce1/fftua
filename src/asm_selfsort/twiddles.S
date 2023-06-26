#define  C            %r8
#define  S            %r9
#define  dC           %r10
#define  dS           %r11
#define  rCos         %ymm0
#define  rSin         %ymm1
#define  dCos         %ymm5
#define  dSin         %ymm6
#define  Cos          %ymm13
#define  Sin          %ymm14
#define  Cos2         %ymm15
#define  Sin2         %ymm11

         .global      twiddle_gen_next
         .global      twiddle_gen_init

         .data
         .align       32

one:     .double      1.0
         .double      1.0
         .double      1.0
         .double      1.0
none:    .double      -1.0
         .double      -1.0
         .double      -1.0
         .double      -1.0
zero:    .double      0.0
         .double      0.0
         .double      0.0
         .double      0.0
half:    .double      0.5
         .double      0.5
         .double      0.5
         .double      0.5

         .text

twiddle_gen_next:
         mov          (%rdx), %rcx
         bsf          8(%rdx), %rax
         shr          %rax
         shl          $5, %rax
         lea          32(%rdx), dC
         lea          (dC, %rax), dS
         lea          (dS, %rax), C
         lea          (C, %rax), S
         inc          %rcx
         mov          %rcx, (%rdx)
         bsf          %rcx, %rax
         shr          %rax
         shl          $5, %rax
         vmovapd      (C, %rax), Cos
         vmovapd      (S, %rax), Sin
         vmovapd      (dC, %rax), dCos
         vmovapd      (dS, %rax), dSin
         vmulpd       dSin, Sin, rCos
         vmulpd       dCos, Sin, rSin
         vfmsub231pd  dCos, Cos, rCos
         vfmadd231pd  dSin, Cos, rSin
         vmovapd      rCos, (%rdi)
         vmovapd      rSin, (%rsi)
L0:      vmovapd      rCos, (C, %rax)
         vmovapd      rSin, (S, %rax)
         sub          $32, %rax
         jns          L0
         ret

twiddle_gen_init:
         xor          %rax, %rax
         mov          %rax, (%rdi)
         mov          %rsi, 8(%rdi)
         bsf          %rsi, %rsi
         shr          %rsi
         shl          $5, %rsi
         lea          32(%rdi), dC
         lea          (dC, %rsi), dS
         lea          (dS, %rsi), C
         lea          (C, %rsi), S
         vmovapd      zero, %ymm0
         vmovapd      one, %ymm1
         vmovapd      half, %ymm12
         vmovapd      none, %ymm10
         vmovapd      %ymm0, Cos
         vmovapd      %ymm10, Sin
         sub          $32, %rsi
L100:    vmovapd      Cos, (dC, %rsi)
         vmovapd      Sin, (dS, %rsi)
         vmovapd      %ymm1, (C, %rsi)
         vmovapd      %ymm0, (S, %rsi)
         sub          $32, %rsi
         js           L200
         vaddpd       %ymm10, Cos, Sin2
         vmulpd       %ymm12, Sin2, Sin2
         vsqrtpd      Sin2, Sin2
         vmovapd      Sin2, Sin
         vaddpd       %ymm1, Cos, Cos2
         vmulpd       %ymm12, Cos2, Cos2
         vsqrtpd      Cos2, Cos2
         vmovapd      Cos2, Cos
         vaddpd       %ymm10, Cos, Sin2
         vmulpd       %ymm12, Sin2, Sin2
         vsqrtpd      Sin2, Sin2
         vmovapd      Sin2, Sin
         vaddpd       %ymm1, Cos, Cos2
         vmulpd       %ymm12, Cos2, Cos2
         vsqrtpd      Cos2, Cos2
         vmovapd      Cos2, Cos
         vmovapd      Cos2, Sin
         vmulpd       Sin, Sin, Sin
         vsubpd       Sin, %ymm1, Sin
         vsqrtpd      Sin, Sin
         vsubpd       Sin, %ymm0, Sin
         jmp          L100
L200:    ret

