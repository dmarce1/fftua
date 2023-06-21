#define    T2            -8(%rbp)
#define    T1            -16(%rbp)
#define    Tptr          -16(%rbp)
#define    C1            -24(%rbp)
#define    S1            -32(%rbp)
#define    C2            -40(%rbp)
#define    S2            -48(%rbp)
#define    N1            -56(%rbp)
#define    N2            -64(%rbp)
#define    N             -72(%rbp)
#define    X             -80(%rbp)
#define    STACK_SIZE    $80

           .global       fft_6step

           .text

fft_6step: enter         STACK_SIZE, $0
           push          %rbx
           mov           %rdi, X
           mov           %rsi, C1
           mov           %rdx, S1
           mov           %rcx, C2
           mov           %r8, S2
           mov           %r9, N1
           mov           16(%rbp), %rax
           mov           %rax, N2
           imul          N1, %rax
           mov           N, %rax
           mov           X, %rdi
           mov           N2, %rsi
           mov           N1, %rdx
           call          fft_scramble_hi
           mov           X, %rdi
           mov           C2, %rsi
           mov           S2, %rdx
           mov           N2, %rcx
           mov           N1, %r8
           call          fft_recursive_hi
           mov           X, %rdi
           mov           C2, %rsi
           mov           S2, %rdx
           mov           N1, %rcx
           mov           N2, %r8
           call          fft_twiddles
           mov           X, %rdi
           mov           C2, %rsi
           mov           S2, %rdx
           mov           N1, %rcx
           mov           N2, %r8
           lea           Tptr, %r9
           call          fft_no2_pre
           mov           N1, %rcx
           imul          N2, %rcx
           bsf           %rcx, %rcx
           mov           %rcx, %rdx
           shr           %rcx
           mov           $1, %rax
           shl           %rcx, %rax
           mov           %rax, %rsi
           mov           %rdx, %rcx
           imul          %rax, %rax
           bsf           %rax, %rax
           sub           %rax, %rcx
           mov           $1, %rbx
           shl           %rcx, %rbx
           mov           %rbx, %rdx
           mov           X, %rdi
           call          fft_transpose_recursive
           mov           X, %rdi
           mov           N2, %rsi
           mov           N1, %rdx
           call          fft_scramble_hi
           mov           X, %rdi
           mov           C1, %rsi
           mov           S1, %rdx
           mov           N2, %rcx
           mov           N1, %r8
           call          fft_recursive_hi
           mov           X, %rdi
           mov           N1, %rsi
           mov           N2, %rdx
           mov           T1, %rcx
           mov           T2, %r8
           call          fft_no2_post
           mov           X, %rdi
           mov           N1, %rsi
           mov           N2, %rdx
           call          fft_real2complex
           mov           X, %rdi
           mov           N1, %rsi
           mov           N2, %rdx
           call          fft_transpose_imag
           pop           %rbx
           leave
           ret
