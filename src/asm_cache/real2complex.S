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
           mov           $1, k1
k1loop:    mov           $1, k2
k2loop:    mov           N2, %rdx
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
           vmovapd       xry, (X, %rcx, 8)
           vmovapd       xiy, (X, %rdx, 8)
           add           $4, k2
           cmp           k2, N2
           jne           k2loop
           inc           k1
           cmp           k1, N1
           jne           k1loop
           pop           %rbx
           pop           %r13
           pop           %r12
           ret

