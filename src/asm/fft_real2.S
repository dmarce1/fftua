#define       X              %r15
#define       N              %r14
#define       No2            %r13

              .global        fft_real2

              .text

fft_real2:    push           %r14
              push           %r15
              mov            %rdi, X
              mov            %rsi, N
              shr            %rsi
              mov            %rsi, No2
              mov            X, %rdi
              mov            No2, %rsi
              call           fft_complex2
              mov            X, %rdi
              mov            N, %rsi
              call           scramble
              mov            X, %rdi
              mov            N, %rsi
              call           complex2real
              pop            %r15
              pop            %r14
              ret
