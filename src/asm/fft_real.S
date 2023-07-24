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


              .align         32
TW45:         .double        0.70710678118654752440
              .double        0.70710678118654752440
              .double        0.70710678118654752440
              .double        0.70710678118654752440
NONE:         .double        -1.0
              .double        -1.0
              .double        -1.0
              .double        -1.0
TWO:          .double        2.0
              .double        2.0
              .double        2.0
              .double        2.0
