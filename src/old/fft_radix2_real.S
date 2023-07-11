#define       X              %r15
#define       N              %r14
#define       No4            %r13
#define       STACK_SIZE     $32
#define       X0             -8(%rbp)
#define       X1             -16(%rbp)
#define       Y0             -24(%rbp)
#define       Y1             -32(%rbp)



              .global        fft_radix2_real


              .text

fft_radix2_real:
              push           %rbp
              mov            %rsp, %rbp
              sub            STACK_SIZE, %rsp
              push           %r15
              push           %r14
              push           %r13
              mov            %rdi, X
              mov            %rsi, N
             // call           scramble
              mov            N, No4
              shr            $2, No4
              lea            (X, No4, 8), %rax
              lea            (%rax, No4, 8), %rcx
              lea            (%rcx, No4, 8), %rdx
              mov            X, X0
              mov            %rax, Y0
              mov            %rcx, X1
              mov            %rdx, Y1
              mov            X0, %rdi
              mov            Y0, %rsi
              mov            No4, %rdx
              call           fft_radix2_complex
              mov            X1, %rdi
              mov            Y1, %rsi
              mov            No4, %rdx
              call           fft_radix2_complex
              mov            X0, %rdi
              mov            Y0, %rsi
              mov            X1, %rdx
              mov            Y1, %rcx
              mov            N, %r8
              call           complex2real
              pop            %r13
              pop            %r14
              pop            %r15
              mov            %rbp, %rsp
              pop            %rbp
              ret
