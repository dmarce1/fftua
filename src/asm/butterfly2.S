#define  SIMD_SIZE      $4
#define  STACK_SIZE     $0
#define  X              %r8
#define  W              %r9
#define  N2             %r10
#define  NLO            %r11
#define  N2o2           %r12
#define  N2NLO          %r13
#define  nlo            %r14
#define  k2             %r15
#define  er0            %ymm0
#define  er1            %ymm1
#define  ei0            %ymm2
#define  ei1            %ymm3
#define  tr0            %ymm4
#define  tr1            %ymm5
#define  ti0            %ymm6
#define  ti1            %ymm7
#define  cos            %ymm8
#define  sin            %ymm9


         .global        fft_swap

fft_butterfly2:
         enter          STACK_SIZE, $0
         push           %r12
         push           %r13
         push           %r14
         push           %r15
         mov            %rdi, X
         mov            %rsi, W
         mov            %rdx, NLO
         mov            %rcx, N2
         mov            N2, N2o2
         shr            N2o2
         mov            N2, %rax
         mul            NLO
         mov            %rax, N2NLO
         mov            $1, k2
k2loop_b:
         cmp            k2, N2o2
         jle            k2loop_e
         mov            k2, %rax
         mul            NLO
         shl            %rax
         vbroadcastsd   (W, %rax, 8), cos
         vbroadcastsd   8(W, %rax, 8), sin
         xor            nlo, nlo
         mov            N2, %rax
         sub            k2, %rax
         sub            k2, %rax
         mul            NLO
         mov            %rax, %rdi
         mov            k2, %rax
         mul            NLO
         mov            %rax, %rcx
nloloop:
         mov            %rcx, %rax
         add            nlo, %rax
         lea            (X, %rax, 8), %rax
         lea            (%rax, N2NLO, 8), %rbx
         vmovapd        (%rax), er0
         vmovapd        (%rbx), er1
         vmovapd        (%rax, %rdi, 8), ei0
         vmovapd        (%rbx, %rdi, 8), ei1
         vmulpd         sin, ei1, tr1
         vmulpd         cos, ei1, ti1
         vfmsub231pd    cos, er1, tr1
         vfmadd231pd    sin, er1, ei1
         vmovapd        tr1, er1
         vaddpd         er1, er0, tr0
         vaddpd         ei1, ei0, ti0
         vsubpd         er1, er0, tr1
         vsubpd         ei1, ei0, ti1
         vmovapd        tr0, (%rax)
         vmovapd        tr1, (%rbx)
         vmovapd        ti1, (%rax, %rdi, 8)
         vmovapd        ti0, (%rbx, %rdi, 8)
         add            SIMD_SIZE, nlo
         cmp            nlo, NLO
         jne            nloloop
         inc            k2
         jmp            k2loop_b
k2loop_e:

         pop            %r15
         pop            %r14
         pop            %r13
         pop            %r12
         leave
         ret
