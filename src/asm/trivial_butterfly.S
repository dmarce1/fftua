#define  SIMD_SIZE    $4
#define  X            %rdi
#define  N            %rsi
#define  M            %rcx
#define  I            %rdx
#define  e0           %ymm0
#define  e1           %ymm1
#define  e2           %ymm2
#define  e3           %ymm3
#define  t0           %ymm4
#define  t1           %ymm5
#define  t2           %ymm6
#define  t3           %ymm7

         .global      trivial_butterfly

         .text

trivial_butterfly:
         mov          N, M
         shr          $2, M
         xor          I, I
loop:
         lea          (X, I, 8), %r8
         lea          (%r8, M, 8), %r9
         lea          (%r9, M, 8), %r10
         lea          (%r10, M, 8), %r11
         vmovapd      (%r8), e0
         vmovapd      (%r9), e1
         vmovapd      (%r10), e2
         vmovapd      (%r11), e3
         vaddpd       e2, e0, t0
         vsubpd       e2, e0, t2
         vaddpd       e3, e1, t1
         vsubpd       e1, e3, t3
         vaddpd       t1, t0, e0
         vsubpd       t1, t0, e2
         vmovapd      e0, (%r8)
         vmovapd      t2, (%r9)
         vmovapd      e2, (%r10)
         vmovapd      t3, (%r11)
         add          SIMD_SIZE, I
         cmp          I, M
         jne          loop
         ret
