#define       X              %r15
#define       N              %r14

              .global        fft_real2


              .text

fft_real2:    push           %r15
              push           %r14
              mov            %rdi, X
              mov            %rsi, N
              shr            N
              mov            X, %rdi
              mov            N, %rsi
              call           fft_core2
              shl            N
              mov            X, %rdi
              mov            N, %rsi
              call           scramble
              mov            X, %rdi
              mov            N, %rsi
              call           complex2real
              pop            %r14
              pop            %r15
              ret
