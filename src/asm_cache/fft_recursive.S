#define  M_SQRT1_2    0.70710678118654752440
#define  N1           $1
#define  N2           %r8
#define  k2           %r9
#define  C            %r10
#define  S            %r11
#define  X            %r12
#define  N2o2         %r13
#define  NTW          %r14
#define  er0          %ymm0
#define  er1          %ymm1
#define  er2          %ymm2
#define  er3          %ymm3
#define  ei0          %ymm4
#define  ei1          %ymm5
#define  ei2          %ymm6
#define  ei3          %ymm7
#define  tr0          %ymm8
#define  tr1          %ymm9
#define  tr2          %ymm10
#define  tr3          %ymm11
#define  ti0          %ymm12
#define  ti1          %ymm13
#define  ti2          %ymm14
#define  ti3          %ymm15
#define  I1           -32(%rbp)
#define  I2           -64(%rbp)
#define  I3           -96(%rbp)
#define  STACK_SIZE   $96

         .global      fft_recursive

         .data

         .align       32
tw45:    .double      M_SQRT1_2
         .double      M_SQRT1_2
         .double      M_SQRT1_2
         .double      M_SQRT1_2
none:    .double      -1.0
         .double      -1.0
         .double      -1.0
         .double      -1.0
iota:    .quad        0
         .quad        1
         .quad        2
         .quad        3
two:     .quad        2
         .quad        2
         .quad        2
         .quad        2
three:   .quad        3
         .quad        3
         .quad        3
         .quad        3

         .text

fft_recursive:
          enter       STACK_SIZE, $0
          push        %rbx
          push        %r12
          push        %r13
          push        %r14
          push        %r15
          mov         %rcx, N2
          mov         %rdx, S
          mov         %rsi, C
          mov         %rdi, X
          shr         $2, N2
          mov         $1, NTW
          call        recurse
          pop         %r15
          pop         %r14
          pop         %r13
          pop         %r12
          pop         %rbx
          leave
          ret

recurse:
          push        X
          shl         $2, NTW
          shr         $2, N2
          call        recurse
          shl         $2, N2
          lea         (X, N2, 8), X
          shr         $2, N2
          call        recurse
          shl         $2, N2
          lea         (X, N2, 8), X
          shr         $2, N2
          call        recurse
          shl         $2, N2
          lea         (X, N2, 8), X
          shr         $2, N2
          call        recurse
          shl         $2, N2
          shr         $2, NTW
          pop         X



done:     ret




