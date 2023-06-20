#define    X             %r8
#define    N1            %r9
#define    N2            %r10
#define    N             %r11
#define    k1            %r12
#define    k2            %r13
#define    N1o2          %r14
#define    N2o2          %r15
#define    prmt_cntrl    $27


           .global       fft_transpose_imag

           .data

           .align        32
none:      .double       -1.0, -1.0, -1.0, -1.0

fft_transpose_imag:
           push          %r12
           push          %r13
           push          %r14
           push          %r15
           push          %rbx
           mov           %rdi, X
           mov           %rsi, N1
           mov           %rdx, N2
           vmovapd       none, %ymm3
           xor           k1, k1
k1loop:    mov           N2o2, k2
           inc           k2
k2loop:    mov           N2, %rax
           imul          k1, %rax
           mov           k1, %rdx
           add           k2, %rax
           add           k2, %rdx
           mov           N, %rbx
           sub           %rdx, %rbx
           mov           %rbx, %rcx
           sub           $3, %rbx
           add           k2, %rcx
           add           k2, %rcx
           vmovapd       (X, %rax, 8), %ymm0
           vmovupd       (X, %rbx, 8), %ymm1
           vmovapd       (X, %rcx, 8), %ymm2
           vmulpd        %ymm3, %ymm2, %ymm2
           vmovapd       %ymm1, (X, %rax, 8)
           vmovupd       %ymm2, (X, %rbx, 8)
           vmovapd       %ymm0, (X, %rcx, 8)
           add           $4, k2
           cmp           N2, k2
           jne           k2loop
           inc           k1
           jne           k1loop
           pop           %rbx
           pop           %r15
           pop           %r14
           pop           %r13
           pop           %r12
           ret


/*for (int k1 = 0; k1 < Mo2; k1++) {
		for (int k2 = Mo2 + 1; k2 < M; k2++) {
			const int i = M * k1 + k2;
			const int j = N - M * k1 - k2;
			const int k = N - M * k1 - (M - k2);
			std::swap(x[i], x[j]);
			std::swap(x[j], x[k]);
			x[j] = -x[j];
		}
	}*/
