#define       X              %r12
#define       N              %r13
#define       C              %r14
#define       N1             %r15
#define       N2             %rbx

              .global        fft_recursive
              .global        none
              .global        two
              .global        tw45

              .data

              .align         32
none:         .double        -1.0
              .double        -1.0
              .double        -1.0
              .double        -1.0
two:          .double        2.0
              .double        2.0
              .double        2.0
              .double        2.0
tw45:         .double        0.70710678118654752440
              .double        0.70710678118654752440
              .double        0.70710678118654752440
              .double        0.70710678118654752440


              .text


fft_recursive:push           %r12
              push           %r13
              push           %r14
              push           %r15
              push           %rbx
              mov            %rdi, X
              mov            %rsi, C
              mov            %rdx, N
              call           dit_nr_entry
              mov            X, %rdi
              mov            N, %rsi
              call           scramble
              mov            X, %rdi
              mov            C, %rsi
              mov            N, %rdx
              call           dit_rn_final
done:         pop            %rbx
              pop            %r15
              pop            %r14
              pop            %r13
              pop            %r12
              ret



























