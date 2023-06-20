#define    C             %r8
#define    S             %r9
#define    X             %r10
#define    N1            %r11
#define    N2            %r12
#define    n1            %r13
#define    k2            %r14
#define    cos0          %ymm0
#define    cos1          %ymm1
#define    cos2          %ymm2
#define    sin0          %ymm3
#define    sin1          %ymm4
#define    sin2          %ymm5
#define    er            %ymm6
#define    ei            %ymm7
#define    tr            %ymm8
#define    ti            %ymm9
#define    prmt_cntrl    $27

           .global       fft_twiddles


           .data

           .align        32
iota:      .quad         0, 1, 2, 3
on:        .double       -1.0, -1.0, -1.0, -1.0

           .text

fft_twiddles:
           push          %r13
           push          %r14
           mov           %r8, N2
           mov           %rcx, N1
           mov           %rdx, S
           mov           %rsi, C
           mov           %rdi, X
           xor           k2, k2
k2loop:    mov           $1, n1
           vmovapd       (C, k2, 8), cos0
           vmovapd       (S, k2, 8), sin0
           vmovapd       cos0, cos1
           vmovapd       sin0, sin1
n1loop:    vmulpd        cos0, sin1, cos2
           vmulpd        sin0, sin1, sin2
           vfmsub231pd   sin0, cos1, cos2
           vfmadd231pd   cos0, cos1, sin2
           vmovapd       cos2, cos1
           vmovapd       sin2, sin1
           mov           n1, %rax
           imul          N2, %rax
           mov           %rax, %rdx
           add           k2, %rax
           add           N2, %rdx
           sub           k2, %rdx
           sub           $3, %rdx
           vmovapd       (X, %rax, 8), er
           vpermpd       prmt_cntrl, (X, %rdx, 8), ei
           vmulpd        cos1, ei, tr
           vmulpd        sin1, ei, ti
           vfmsub231pd   sin1, er, tr
           vfmadd231pd   cos1, er, ti
           vmovapd       tr, er
           vpermpd       prmt_cntrl, ti, ei
           vmovapd       er, (X, %rax, 8)
           vmovupd       ei, (X, %rdx, 8)
           inc           n1
           cmp           n1, N1
           jne           n1loop
           add           $4, k2
           cmp           k2, N2
           jne           k2loop
           pop           %r14
           pop           %r13
           ret

