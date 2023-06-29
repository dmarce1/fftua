#define       STACK_SIZE     $512
#define       ER0            -32(%rbp)
#define       ER1            -64(%rbp)
#define       ER2            -96(%rbp)
#define       ER3            -128(%rbp)
#define       EI0            -160(%rbp)
#define       EI1            -192(%rbp)
#define       EI2            -224(%rbp)
#define       EI3            -256(%rbp)
#define       OR0            -288(%rbp)
#define       OR1            -320(%rbp)
#define       OR2            -352(%rbp)
#define       OR3            -384(%rbp)
#define       OI0            -416(%rbp)
#define       OI1            -448(%rbp)
#define       OI2            -480(%rbp)
#define       OI3            -512(%rbp)
#define       er0            %ymm0
#define       er1            %ymm1
#define       er2            %ymm2
#define       er3            %ymm3
#define       ei0            %ymm4
#define       ei1            %ymm5
#define       ei2            %ymm6
#define       ei3            %ymm7
#define       or0            %ymm8
#define       or1            %ymm9
#define       or2            %ymm10
#define       or3            %ymm11
#define       oi0            %ymm12
#define       oi1            %ymm13
#define       oi2            %ymm14
#define       oi3            %ymm15


              .global        butterfly8
              .global        butterfly8_test

              .data

              .align         32
W0:           .double        0.70710678118654752440
              .double        0.70710678118654752440
              .double        0.70710678118654752440
              .double        0.70710678118654752440
W1:           .double        -0.70710678118654752440
              .double        -0.70710678118654752440
              .double        -0.70710678118654752440
              .double        -0.70710678118654752440


              .text

butterfly8:   push           %rbp
              mov            %rsp, %rbp
              sub            STACK_SIZE, %rsp
              vmovupd        er0, ER0
              vmovupd        er1, ER1
              vmovupd        er2, ER2
              vmovupd        er3, ER3
              vmovupd        ei0, EI0
              vmovupd        ei1, EI1
              vmovupd        ei2, EI2
              vmovupd        ei3, EI3
              vaddpd         or2, or0, er0
              vaddpd         oi2, oi0, ei0
              vsubpd         or2, or0, er2
              vsubpd         oi2, oi0, ei2
              vaddpd         or3, or1, er1
              vaddpd         oi3, oi1, ei1
              vsubpd         or3, or1, er3
              vsubpd         oi3, oi1, ei3
              vaddpd         er1, er0, or0
              vaddpd         ei1, ei0, oi0
              vaddpd         ei3, er2, or1
              vsubpd         er3, ei2, oi1
              vsubpd         er1, er0, or2
              vsubpd         ei1, ei0, oi2
              vsubpd         ei3, er2, or3
              vaddpd         er3, ei2, oi3
              vmovupd        or0, OR0
              vmovupd        or1, OR1
              vmovupd        or2, OR2
              vmovupd        or3, OR3
              vmovupd        oi0, OI0
              vmovupd        oi1, OI1
              vmovupd        oi2, OI2
              vmovupd        oi3, OI3
              vmovupd        ER0, er0
              vmovupd        ER1, er1
              vmovupd        EI0, ei0
              vmovupd        EI1, ei1
              vaddpd         ER2, er0, or0
              vaddpd         EI2, ei0, oi0
              vsubpd         ER2, er0, or2
              vsubpd         EI2, ei0, oi2
              vaddpd         ER3, er1, or1
              vaddpd         EI3, ei1, oi1
              vsubpd         ER3, er1, or3
              vsubpd         EI3, ei1, oi3
              vaddpd         or1, or0, er0
              vaddpd         oi1, oi0, ei0
              vaddpd         oi3, or2, er1
              vsubpd         or3, oi2, ei1
              vsubpd         or1, or0, er2
              vsubpd         oi1, oi0, ei2
              vsubpd         oi3, or2, er3
              vaddpd         or3, oi2, ei3
              vmovupd        OR1, or1
              vmovupd        OI3, or3
              vmovupd        OR1, oi1
              vmovupd        OR3, oi3
              vaddpd         OI1, or1, or1
              vsubpd         OR3, or3, or3
              vsubpd         OI1, oi1, oi1
              vaddpd         OI3, oi3, oi3
              vmulpd         W0, or1, or1
              vmulpd         W0, or3, or3
              vmulpd         W1, oi1, oi1
              vmulpd         W1, oi3, oi3
              vmovupd        or1, OR1
              vmovupd        or3, OR3
              vmovupd        oi1, OI1
              vmovupd        oi3, OI3
              vsubpd         OR0, er0, or0
              vsubpd         OI0, ei0, oi0
              vaddpd         OR0, er0, er0
              vaddpd         OI0, ei0, ei0
              vsubpd         OR1, er1, or1
              vsubpd         OI1, ei1, oi1
              vaddpd         OR1, er1, er1
              vaddpd         OI1, ei1, ei1
              vsubpd         OI2, er2, or2
              vaddpd         OR2, ei2, oi2
              vaddpd         OI2, er2, er2
              vsubpd         OR2, ei2, ei2
              vsubpd         OR3, er3, or3
              vsubpd         OI3, ei3, oi3
              vaddpd         OR3, er3, er3
              vaddpd         OI3, ei3, ei3
              mov            %rbp, %rsp
              pop            %rbp
              ret


butterfly8_test:
              vbroadcastsd   0(%rdi), er0
              vbroadcastsd   8(%rdi), ei0
              vbroadcastsd   16(%rdi), or0
              vbroadcastsd   24(%rdi), oi0
              vbroadcastsd   32(%rdi), er1
              vbroadcastsd   40(%rdi), ei1
              vbroadcastsd   48(%rdi), or1
              vbroadcastsd   56(%rdi), oi1
              vbroadcastsd   64(%rdi), er2
              vbroadcastsd   72(%rdi), ei2
              vbroadcastsd   80(%rdi), or2
              vbroadcastsd   88(%rdi), oi2
              vbroadcastsd   96(%rdi), er3
              vbroadcastsd   104(%rdi), ei3
              vbroadcastsd   112(%rdi), or3
              vbroadcastsd   120(%rdi), oi3
              call           butterfly8
              vmovq          %xmm0, 0(%rdi)
              vmovq          %xmm4, 8(%rdi)
              vmovq          %xmm1, 16(%rdi)
              vmovq          %xmm5, 24(%rdi)
              vmovq          %xmm2, 32(%rdi)
              vmovq          %xmm6, 40(%rdi)
              vmovq          %xmm3, 48(%rdi)
              vmovq          %xmm7, 56(%rdi)
              vmovq          %xmm8, 64(%rdi)
              vmovq          %xmm12, 72(%rdi)
              vmovq          %xmm9, 80(%rdi)
              vmovq          %xmm13, 88(%rdi)
              vmovq          %xmm10, 96(%rdi)
              vmovq          %xmm14, 104(%rdi)
              vmovq          %xmm11, 112(%rdi)
              vmovq          %xmm15, 120(%rdi)
              ret

