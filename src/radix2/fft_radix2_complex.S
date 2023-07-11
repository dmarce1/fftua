#define       X              %r15
#define       stack_size     %r14
#define       N              %r13
#define       N2             %r12
#define       W              %r11
#define       k2             %r9
#define       Wptr           %r8
#define       e0             %ymm0
#define       e1             %ymm1
#define       e2             %ymm2
#define       e3             %ymm3
#define       t0             %ymm4
#define       t1             %ymm5
#define       t2             %ymm6
#define       t3             %ymm7
#define       W1             %ymm8
#define       W2             %ymm9
#define       W3             %ymm10
#define       tmp0           %ymm11
#define       tmp1           %ymm12
#define       tmp2           %ymm13
#define       tmp3           %ymm14
#define       rot90          %ymm15


              .text

fft_radix2_complex:
              push           %rbx
              push           %r15
              push           %r14
              push           %r13
              push           %r12
              mov            %rdi, X
              mov            %rsi, N
              bsf            N, %rax
              inc            %rax
              shl            $3, %rax
              push           %rbp
              mov            %rsp, %rbp
              sub            %rax, %rsp
              mov            %rsp, Wptr
              mov            %rax, stack_size
              bsf            N, %rdx
              mov            $2, %rdi
twiddle_loop: mov            %rdi, %rcx
              push           %rdi
              mov            $1, %rdi
              shl            %rcx, %rdi
              push           %rdx
              push           %r8
              call           twiddles_fwd_complex2
              pop            %r8
              pop            %rcx
              pop            %rdi
              mov            %rax, (Wptr, %rdi, 8)
              mov            %rcx, %rdx
              cmp            %rdi, %rdx
              je             twiddle_esc
              inc            %rdi
              jmp            twiddle_loop
twiddle_esc:  mov            N, N2
              shr            $2, N2
              call           next_level
              add            stack_size, %rsp
              mov            %rsp, %rbp
              pop            %rbp
              pop            %r12
              pop            %r13
              pop            %r14
              pop            %r15
              pop            %rbx
              ret
next_level:   imul           $2, N2, %rcx
              imul           $3, N2, %rdx
              bsf            N, %rax
              mov            (Wptr, %rax, 8), W
              vmovapd        ROT90, rot90
              xor            k2, k2
k2_loop:      vmovapd        (W, k2, 8), W1
              lea            (X, k2, 8), %rax
              vmovapd        (%rax), e0
              vmovapd        (%rax, N2, 8), e1
              vmovapd        (%rax, %rcx, 8), e2
              vmovapd        (%rax, %rdx, 8), e3
              vaddpd         e2, e0, t0
              vsubpd         e2, e0, t2
              vaddpd         e3, e1, t1
              vsubpd         e3, e1, t3
              vpermpd        $177, t3, t3
              vmulpd         rot90, t3, t3
              vaddpd         t1, t0, e0
              vsubpd         t1, t0, e2
              vaddpd         t3, t2, e1
              vsubpd         t3, t2, e3
              vmulpd         W1, W1, tmp0
              vpermpd        $5, W1, tmp1
              vmulpd         rot90, tmp1, tmp1
              vmulpd         tmp1, W1, tmp2
              vhsubpd        tmp1, tmp2, W2
              vmulpd         W1, W2, tmp0
              vpermpd        $5, W2, tmp1
              vmulpd         rot90, tmp1, tmp1
              vmulpd         tmp1, W1, tmp2
              vhsubpd        tmp1, tmp2, W3
              vmulpd         e1, W1, tmp0
              vpermpd        $5, W1, tmp1
              vmulpd         rot90, tmp1, tmp1
              vmulpd         tmp1, e1, tmp2
              vhsubpd        tmp1, tmp2, e1
              vmulpd         e2, W2, tmp0
              vpermpd        $5, W2, tmp1
              vmulpd         rot90, tmp1, tmp1
              vmulpd         tmp1, e2, tmp2
              vhsubpd        tmp1, tmp2, e2
              vmulpd         e3, W3, tmp0
              vpermpd        $5, W3, tmp1
              vmulpd         rot90, tmp1, tmp1
              vmulpd         tmp1, e3, tmp2
              vhsubpd        tmp1, tmp2, e3
              vmovapd        e0, (%rax)
              vmovapd        e2, (%rax, N2, 8)
              vmovapd        e1, (%rax, %rcx, 8)
              vmovapd        e3, (%rax, %rdx, 8)
              add            $4, k2
              cmp            k2, N2
              jne            k2_loop
              cmp            $64, N
              jl             done
              push           X
              imul           $3, N2, %rdx
              lea            (X, %rdx, 8), %rax
              push           %rax
              imul           $2, N2, %rdx
              lea            (X, %rdx, 8), %rax
              push           %rax
              lea            (X, N2, 8), %rax
              push           %rax
              shr            $2, N2
              shr            $2, N
              call           next_level
              pop            X
              call           next_level
              pop            X
              call           next_level
              pop            X
              call           next_level
              pop            X
              shl            $2, N2
              shl            $2, N
done:         ret
              .align         32
ROT90:        .double        1.0
              .double        -1.0
              .double        1.0
              .double        -1.0
TWO:          .double        2.0
              .double        2.0
              .double        2.0
              .double        2.0
tw16_1r:      .double        +1.000000000000000000
              .double        +0.923879532511286756
              .double        +0.707106781186547524
              .double        +0.382683432365089771
tw16_1i:      .double        -0.000000000000000000
              .double        -0.382683432365089771
              .double        -0.707106781186547524
              .double        -0.923879532511286756
tw16_2r:      .double        +1.000000000000000000
              .double        +0.707106781186547524
              .double        +0.000000000000000000
              .double        -0.707106781186547524
tw16_2i:      .double        +0.000000000000000000
              .double        -0.707106781186547524
              .double        -1.000000000000000000
              .double        -0.707106781186547524
tw32_1r:      .double        1.000000000000000000
              .double        0.980785280403230000
              .double        0.923879532511287000
              .double        0.831469612302545000
tw32_2r:      .double        0.707106781186548000
              .double        0.555570233019602000
              .double        0.382683432365090000
              .double        0.195090322016128000
tw32_3r:      .double        0.000000000000000000
              .double        -0.195090322016128000
              .double        -0.382683432365090000
              .double        -0.555570233019602000
tw32_4r:      .double        -0.707106781186547000
              .double        -0.831469612302545000
              .double        -0.923879532511287000
              .double        -0.980785280403230000
tw32_1i:      .double        0.000000000000000000
              .double        -0.195090322016128000
              .double        -0.382683432365090000
              .double        -0.555570233019602000
tw32_2i:      .double        -0.707106781186547000
              .double        -0.831469612302545000
              .double        -0.923879532511287000
              .double        -0.980785280403230000
tw32_3i:      .double        -1.000000000000000000
              .double        -0.980785280403230000
              .double        -0.923879532511287000
              .double        -0.831469612302545000
tw32_4i:      .double        -0.707106781186548000
              .double        -0.555570233019602000
              .double        -0.382683432365090000
              .double        -0.195090322016129000

