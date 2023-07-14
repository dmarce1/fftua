#define       X              %r15
#define       N              %r14

              .global        fft_complex2


              .text

fft_complex2: push           %r15
              push           %r14
              mov            %rdi, X
              mov            %rsi, N
              call           fft_core2
              mov            X, %rdi
              mov            N, %rsi
              call           scramble_cmplx
              pop            %r14
              pop            %r15
              ret
