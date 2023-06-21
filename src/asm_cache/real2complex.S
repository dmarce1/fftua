#define    rptr          %r8
#define    X             %r8
#define    N1            %r9
#define    N2            %r10
#define    k1            %r11
#define    k2            %r12
#define    N             %r13
#define    rex           %ymm0
#define    rey           %ymm1
#define    imx           %ymm2
#define    imy           %ymm3
#define    xrx           %ymm4
#define    xix           %ymm5
#define    xry           %ymm6
#define    xiy           %ymm7
#define    srex           %xmm0
#define    srey           %xmm1
#define    simx           %xmm2
#define    simy           %xmm3
#define    sxrx           %xmm4
#define    sxix           %xmm5
#define    sxry           %xmm6
#define    sxiy           %xmm7
#define    prmt_cntrl    $27


           .global       fft_real2complex


           .text

fft_real2complex:
           push          %r12
           push          %r13
           push          %rbx
           mov           %rdi, X
           mov           %rsi, N1
           mov           %rdx, N2
           mov           N1, N
           imul          N2, N
           mov           $1, k1
k1loop:    mov           $1, k2
scalar:    mov           N2, %rdx
           imul          k1, %rdx
           mov           k2, %rax
           add           %rdx, %rax
           mov           N, %rbx
           sub           %rax, %rbx
           mov           %rax, %rcx
           sub           k2, %rcx
           sub           k2, %rcx
           mov           N, %rdx
           sub           %rcx, %rdx
           movsd         (X, %rax, 8), srex
           movsd         (X, %rbx, 8), srey
           movsd         (X, %rcx, 8), simx
           movsd         (X, %rdx, 8), simy
           movsd         srex, sxrx
           movsd         srey, sxix
           movsd         srex, sxry
           movsd         simx, sxiy
           vsubsd        simy, sxrx, sxrx
           vaddsd        simx, sxix, sxix
           vaddsd        simy, sxry, sxry
           vsubsd        srey, sxiy, sxiy
           movsd         sxrx, (X, %rax, 8)
           movsd         sxix, (X, %rbx, 8)
           movsd         sxry, (X, %rcx, 8)
           movsd         sxiy, (X, %rdx, 8)
           inc           k2
           cmp           $4, k2
           jne           scalar
simd_b:    cmp           k2, N2
           jle           simd_e
		   mov           N2, %rdx
           imul          k1, %rdx
           mov           k2, %rax
           add           %rdx, %rax
           mov           N, %rbx
           sub           %rax, %rbx
           mov           %rax, %rcx
           sub           k2, %rcx
           sub           k2, %rcx
           mov           N, %rdx
           sub           %rcx, %rdx
           sub           $3, %rcx
           sub           $3, %rdx
           vmovapd       (X, %rax, 8), rex
           vmovapd       (X, %rbx, 8), rey
           vpermpd       prmt_cntrl, (X, %rcx, 8), imx
           vpermpd       prmt_cntrl, (X, %rdx, 8), imy
           vmovapd       rex, xrx
           vmovapd       rey, xix
           vmovapd       rex, xry
           vmovapd       imx, xiy
           vsubpd        imy, xrx, xrx
           vaddpd        imx, xix, xix
           vaddpd        imy, xry, xry
           vsubpd        rey, xiy, xiy
           vpermpd       prmt_cntrl, xry, xry
           vpermpd       prmt_cntrl, xiy, xiy
           vmovapd       xrx, (X, %rax, 8)
           vmovapd       xix, (X, %rbx, 8)
           vmovupd       xry, (X, %rcx, 8)
           vmovupd       xiy, (X, %rdx, 8)
           add           $4, k2
           jmp           simd_b
simd_e:    inc           k1
           cmp           k1, N1
           jne           k1loop
           pop           %rbx
           pop           %r13
           pop           %r12
           ret

