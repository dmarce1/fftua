#define SIMD_SIZE      $4
#define SIMD_SIZE2     $8
#define X              %rdi
#define W              %rsi
#define N              %r10
#define k2             %r11
#define N2             %r12
#define n1             %r9
#define No2            %r8
#define er             %ymm0
#define ei             %ymm1
#define ti             %ymm2
#define tr             %ymm3

       .global         fft_twiddles

       .text
fft_twiddles:
       push            %r12
       movq            %rdx, N
       movq            N, No2
       shr             No2
       movq            N, %rax
       mul             %rax
       movq            %rax, N2
       mov             $1, k2
oloop:
       mov             N2, %rcx
       mov             k2, %rax
       mul             N
       sub             %rax, %rcx
       mov             %rax, %rdx
       shl             %rax
       xor             n1, n1
iloop:
       vmovapd         (X, %rcx, 8), ei
       vmovapd         (X, %rdx, 8), er
       vmulpd          32(W, %rax, 8), ei, tr
       vmulpd          (W, %rax, 8), ei, ti
       vfmsub231pd     (W, %rax, 8), er, tr
       vfmadd231pd     32(W, %rax, 8), er, ti
       vmovapd         tr, (X, %rdx, 8)
       vmovapd         ti, (X, %rcx, 8)
       add             SIMD_SIZE, n1
       add             SIMD_SIZE, %rdx
       add             SIMD_SIZE, %rcx
       add             SIMD_SIZE2, %rax
       cmp             N, n1
       jl              iloop
       inc             k2
       cmp             No2, k2
       jl              oloop
       pop             %r12
       emms
       ret

