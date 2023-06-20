// void fft_no2_post( double* X, int N1, int N2, double T1, double T2 );
// void fft_no2_pre( double* X, const double* C, const double* S, int N1, int N2, double* rc );
// void fft_recursive_hi( double* X, const double* C, const double* S, int N, int M );
// void fft_transpose_imag( double* X, int N1, int N2 );
// void fft_real2complex( double* X, int N1, int N2 );
// void fft_transpose( double* X, const double* C, const double* S, int N1, int N2 );
// void fft_twiddles( double* X, const double* C, const double* S, int N1, int N2 );

#define    X             %r10
#define    C             %r11
#define    S             %r12
#define    N             %r13
#define    N1            %r14
#define    N2            %r15
#define    T1            -8(%rbp)
#define    T2            -16(%rbp)
#define    Tptr          -16(%rbp)
#define    STACK_SIZE    $16

           .global       fft_6step

           .text

fft_6step: push          %r12
           push          %r13
           push          %r14
           push          %r15
           mov           %rdi, X
           mov           %rsi, N
           mov           %rdx, C
           mov           %rcx, S
           bsf           N, %rax
           mov           %rax, %rbx
           shr           %rbx
           sub           %rbx, %rax
           mov           $1, N1
           mov           %rax, %rcx
           shl           %rcx, N1
           mov           $1, N2
           mov           %rax, %rcx
           shl           %rcx, N2
           mov           X, %rdi
           mov           N2, %rsi
           mov           N1, %rdx
           call          fft_scramble_hi
           mov           X, %rdi
           mov           C, %rsi
           mov           S, %rdx
           mov           N2, %rcx
           mov           N1, %r8
           call          fft_recursive_hi
           mov           X, %rdi
           mov           C, %rsi
           mov           S, %rdx
           mov           N1, %rcx
           mov           N2, %r8
           call          fft_twiddles
           mov           X, %rdi
           mov           C, %rsi
           mov           S, %rdx
           mov           N1, %rcx
           mov           N2, %r8
           lea           Tptr, %r9
           call          fft_no2_pre
           mov           X, %rdi
           mov           C, %rsi
           mov           S, %rdx
           mov           N1, %rcx
           mov           N2, %r8
           call          fft_transpose_recursive
           mov           X, %rdi
           mov           N2, %rsi
           mov           N1, %rdx
           call          fft_scramble_hi
           mov           X, %rdi
           mov           C, %rsi
           mov           S, %rdx
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
           pop           %r15
           pop           %r14
           pop           %r13
           pop           %r12
           ret
