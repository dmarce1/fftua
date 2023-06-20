#define    rptr          %r8
#define    X             %r9
#define    C             %r10
#define    S             %r11
#define    N1            %r12
#define    N2            %r13
#define    T1            %ymm14
#define    T2            %ymm15
#define    tmp3          -8(%rbp)
#define    tmp2          -16(%rbp)
#define    tmp1          -24(%rbp)
#define    tmp0          -32(%rbp)
#define    STACK_SIZE    $32

           .global       fft_no2_pre

           .data

           .align        32
ntwo:      .double       -2.0, -2.0, -2.0, -2.0
on:        .double       -1.0, -1.0, -1.0, -1.0
evens:     .quad         0, 2, 4, 6
odds:      .quad         1, 3, 5, 7
iota:      .quad         0, 1, 2, 3

           .text

fft_no2_pre:
           enter         STACK_SIZE, $0
           push          %r12
           push          %r13
           mov           %r9, rptr
           mov           %r8, N2
           mov           %rcx, N1
           mov           %rdx, C
           mov           %rsi, S
           mov           %rdi, X
           vpxor         T1, T1, T1
           vpxor         T2, T2, T2
           mov           N2, %rax
           shr           %rax
           imul          N1, %rax
		   lea           (X, %rax, 8), X
		   vmovapd       on, %ymm0
           vmovapd       evens, %ymm1
           vmovapd       odds, %ymm2
           vmovapd       ntwo, %ymm13
           vmovapd       iota, %ymm12
		   xor           %rdi, %rdi
loop0:     lea           (C, %rdi, 8), %rax
           lea           (X, %rdi, 8), %rdx
           vmovapd       %ymm0, %ymm3
           vmovapd       %ymm0, %ymm4
           vmovapd       %ymm0, %ymm5
           vmovapd       %ymm0, %ymm6
           vgatherqpd    %ymm3, (%rax, %ymm1, 8), %ymm7
           vgatherqpd    %ymm4, (%rax, %ymm2, 8), %ymm8
           vgatherqpd    %ymm5, (%rdx, %ymm1, 8), %ymm9
           vgatherqpd    %ymm6, (%rdx, %ymm2, 8), %ymm10
           vfmadd231pd   %ymm9, %ymm7, T1
           vfmadd231pd   %ymm10, %ymm8, T2
	       add           $8, %rdi
	       cmp           N1, %rdi
	       jne           loop0
           vmovupd       T1, tmp0
           movsd         tmp0, %xmm0
           movsd         tmp1, %xmm1
           movsd         tmp2, %xmm2
           movsd         tmp3, %xmm3
           vaddsd        %xmm1, %xmm0, %xmm0
           vaddsd        %xmm2, %xmm0, %xmm0
           vaddsd        %xmm3, %xmm0, %xmm0
           movsd         %xmm0, (rptr)
           vmovupd       T2, tmp0
           movsd         tmp0, %xmm0
           movsd         tmp1, %xmm1
           movsd         tmp2, %xmm2
           movsd         tmp3, %xmm3
           vaddsd        %xmm1, %xmm0, %xmm0
           vaddsd        %xmm2, %xmm0, %xmm0
           vaddsd        %xmm3, %xmm0, %xmm0
           movsd         %xmm0, 8(rptr)
 	       xor           %rdi, %rdi
loop1:     vmovapd       %ymm0, %ymm3
           vmovapd       %ymm0, %ymm4
           lea           (X, %rdi, 8), %rax
           lea           (S, %rdi, 8), %rdx
           vgatherqpd    %ymm3, (%rax, %ymm12, 8), %ymm10
           vgatherqpd    %ymm4, (%rdx, %ymm12, 8), %ymm11
           vmulpd        %ymm11, %ymm10, %ymm10
           vmulpd        %ymm13, %ymm10, %ymm10
           vmovapd       %ymm10, (X, %rdi, 8)
	       add           $4, %rdi
	       cmp           N1, %rdi
	       jne           loop1
           pop           %r13
           pop           %r12
           leave
           ret
