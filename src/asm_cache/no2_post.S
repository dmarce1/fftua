#define    X             %r8
#define    N1            %r9
#define    N2            %r10
#define    k1            %r11
#define    N1o2          %r12
#define    T1            %xmm12
#define    T2            %xmm13
#define    T3            %xmm14
#define    T4            %xmm15


           .global       fft_no2_post

           .data

           .align        8
half:      .double       0.5
nhalf:     .double       -0.5

           .text

fft_no2_post:
           push          %r12
           push          %rbx
           movq          %r8, T2
           movq          %rcx, T1
           mov           %rdx, N2
           mov           %rsi, N1
           mov           %rdi, X
           mov           N1, N1o2
           shr           N1o2
           mov           N2, %rax
           dec           %rax
           imul          N2, %rax
           add           N1o2, %rax
           movq          (X, %rax, 8), T3
           movsd         (X, N1o2, 8), %xmm0
           vmulsd        nhalf, %xmm0, %xmm0
           movsd         %xmm0, (X, %rax, 8)
           vaddsd        T1, T2, %xmm0
           movsd         %xmm0, (X, N1o2, 8)
           mov           $1, k1
loop0:     mov           N1, %rax
           mov           N2, %rcx
           sub           k1, %rax
           mov           %rax, %rbx
           dec           %rax
           imul          k1, %rcx
           imul          N2, %rax
           imul          N2, %rbx
           add           N1o2, %rax
           add           N1o2, %rbx
           add           N1o2, %rcx
           movsd         (X, %rax, 8), T4
           movsd         (X, %rbx, 8), %xmm0
           movsd         (X, %rcx, 8), %xmm1
           movsd         %xmm0, %xmm2
           vsubsd        %xmm1, %xmm0, %xmm0
           movsd         %xmm0, (X, %rax, 8)
           vaddsd        T3, %xmm2, %xmm0
           movsd         %xmm0, (X, %rcx, 8)
           inc           k1
           mov           N1o2, %rax
           dec           %rax
           cmp           k1, %rax
           jne           loop0
           mov           N2, %rdx
           shr           %rbx
           mov           N1o2, %rax
           add           %rbx, %rax
           sub           N1o2, %rbx
           vsubsd        T1, T2, %xmm0
           movsd         %xmm0, (X, %rbx, 8)
           movsd         (X, %rax, 8), %xmm0
           vmulsd        half, %xmm0, %xmm0
           movsd         %xmm0, (X, %rax, 8)
           pop           %rbx
           pop           %r12
           ret

    /*
    T3 = x[Mo2 + M * (M - 1)];
	x[Mo2 + M * (M - 1)] = -0.5 * x[Mo2];
	x[Mo2] = (T1 + T2);
	for (int k1 = 1; k1 < Mo2 - 1; k1++) {
		T4 = x[Mo2 + M * (M - k1 - 1)];
		x[Mo2 + M * (M - k1 - 1)] = x[Mo2 + M * (M - k1)] - x[Mo2 + M * k1];
		x[Mo2 + M * k1] = T3 + x[Mo2 + M * (k1 - 1)];
		T3 = T4;
	}
	x[Mo2 + No2] = 0.5 * x[No2 + Mo2];
	x[No2 - Mo2] = (T1 - T2);*/
