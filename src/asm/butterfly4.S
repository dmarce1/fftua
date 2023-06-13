#define  STACK_SIZE   $112
#define  M_SQRT1_2    0.70710678118654752440
#define  N1           $4
#define  X            %r8
#define  W            %r9
#define  ihi          %r10
#define  N2           %r11
#define  khi          %r12
#define  klo          %r13
#define  k2           %r14
#define  er0          %ymm0
#define  ei0          %ymm1
#define  er1          %ymm2
#define  ei1          %ymm3
#define  er2          %ymm4
#define  ei2          %ymm5
#define  er3          %ymm6
#define  ei3          %ymm7
#define  tr0          %ymm8
#define  tr1          %ymm9
#define  tr2          %ymm10
#define  tr3          %ymm11
#define  ti0          %ymm12
#define  ti1          %ymm13
#define  ti2          %ymm14
#define  ti3          %ymm15
#define  NHI          -8(%rbp)
#define  klo1         -16(%rbp)
#define  N            -24(%rbp)
#define  NHIoN1       -32(%rbp)
#define  N2oN1        -40(%rbp)
#define  N2N1         -48(%rbp)
#define  N2o2         -56(%rbp)
#define  khi1         -64(%rbp)
#define  N2oN1o2      -72(%rbp)
#define  TW1          -80(%rbp)
#define  TW2          -88(%rbp)
#define  TW3          -96(%rbp)
#define  N2o2hi       -104(%rbp)
#define  N2o2lo       -112(%rbp)

         .global      butterfly4

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

         .text

butterfly4:
         enter        STACK_SIZE, $0
         push         %rbx
         push         %r12
         push         %r13
         push         %r14
         mov          %rdi, X
         mov          %rsi, W
         mov          %rdx, NHI
         mov          %rcx, N2
         shr          $2, %rcx
         shr          $2, %rdx
         mov          %rcx, N2oN1
         mov          %rdx, NHIoN1
         mov          N2, %rdx
         mov          NHI, %rax
         mul          %rdx
         shl          $2, %rax
         mov          %rax, N
         mov          N2, %rbx
         shl          $2, %rbx
         mov          %rbx, N2N1
         mov          N2, %rcx
         shr          %rcx
         mov          %rcx, N2o2
         mov          %rcx, %rax
         mov          %rcx, %rbx
         and          $3, %rbx
         shr          $2, %rax
         mov          %rax, N2o2hi
         mov          %rbx, N2o2lo
         shr          $2, %rcx
         cmp          $1, %rcx
         mov          $1, %rax
         cmovl        %rax, %rcx
         mov          %rcx, N2oN1o2
         mov          NHI, %rax
         shl          %rax
         mov          %rax, %rbx
         mov          %rax, %rcx
         add          %rax, %rbx
         add          %rbx, %rcx
         mov          %rax, TW1
         mov          %rbx, TW2
         mov          %rcx, TW3
         xor          ihi, ihi
ihiloop:
         mov          ihi, %rax
         imul         N2N1, %rax
         lea          (X, %rax, 8), %rax
         lea          (%rax, N2, 8), %rbx
         lea          (%rbx, N2, 8), %rcx
         lea          (%rcx, N2, 8), %rdx
         vmovapd      (%rax), er0
         vmovapd      (%rbx), er1
         vmovapd      (%rcx), er2
         vmovapd      (%rdx), er3
         vaddpd       er2, er0, tr0
         vsubpd       er2, er0, tr2
         vaddpd       er3, er1, tr1
         vsubpd       er1, er3, tr3
         vaddpd       tr1, tr0, er0
         vsubpd       tr1, tr0, er2
         vmovapd      er0, (%rax)
         vmovapd      tr2, (%rbx)
         vmovapd      er2, (%rcx)
         vmovapd      tr3, (%rdx)
         cmp          $1, N2
         jle          done
         mov          N2o2lo, %rax
         imul         NHIoN1, %rax
         add          ihi, %rax
         imul         N2, %rax
         add          N2o2hi, %rax
         shl          $2, %rax
         lea          (X, %rax, 8), %rax
         lea          (%rax, N2, 8), %rbx
         lea          (%rbx, N2, 8), %rcx
         lea          (%rcx, N2, 8), %rdx
         vmovapd      (%rax), er0
         vmovapd      (%rbx), er1
         vmovapd      (%rcx), er2
         vmovapd      (%rdx), er3
         vaddpd       er3, er1, tr0
         vsubpd       er3, er1, tr2
         vmulpd       tw45, tr2, tr1
         vmulpd       tw45, tr0, tr3
         vmovapd      er0, tr0
         vmovapd      er2, tr2
         vaddpd       tr1, tr0, er0
         vaddpd       tr3, tr2, er3
         vsubpd       tr1, tr0, er1
         vsubpd       tr3, tr2, er2
         vmulpd       none, er3, er3
         vmovapd      er0, (%rax)
         vmovapd      er1, (%rbx)
         vmovapd      er2, (%rcx)
         vmovapd      er3, (%rdx)
         cmp          $2, N2
         jle          done
         xor          klo, klo
kloloop:
         xor          khi, khi
khiloop:
         mov          khi, k2
         shl          $2, k2
         or           klo, k2
         cmp          k2, N2o2
         jle          khiloop_cont
         test         k2, k2
         je           khiloop_cont
         mov          k2, %rdx
         mulx         TW1, %rax, %rdi
         mulx         TW2, %rbx, %rdi
         mulx         TW3, %rcx, %rdi
         vbroadcastsd 0(W, %rax, 8), tr1
         vbroadcastsd 8(W, %rax, 8), ti1
         vbroadcastsd 0(W, %rbx, 8), tr2
         vbroadcastsd 8(W, %rbx, 8), ti2
         vbroadcastsd 0(W, %rcx, 8), tr3
         vbroadcastsd 8(W, %rcx, 8), ti3
         mov          N2, %rax
         sub          k2, %rax
         mov          %rax, %rdi
         shr          $2, %rax
         and          $3, %rdi
         mov          klo, %rsi
         imul         NHIoN1, %rdi
         imul         NHIoN1, %rsi
         add          ihi, %rdi
         add          ihi, %rsi
         imul         N2, %rdi
         imul         N2, %rsi
         add          %rax, %rdi
         add          khi, %rsi
         shl          $2, %rdi
         shl          $2, %rsi
         sub          %rsi, %rdi
         lea          (X, %rsi, 8), %rax
         lea          (%rax, N2, 8), %rbx
         lea          (%rbx, N2, 8), %rcx
         lea          (%rcx, N2, 8), %rdx
         vmovapd      (%rax), er0
         vmovapd      (%rbx), er1
         vmovapd      (%rcx), er2
         vmovapd      (%rdx), er3
         vmovapd      (%rax, %rdi, 8), ei0
         vmovapd      (%rbx, %rdi, 8), ei1
         vmovapd      (%rcx, %rdi, 8), ei2
         vmovapd      (%rdx, %rdi, 8), ei3
         vmulpd       ti1, ei1, tr0
         vmulpd       tr1, ei1, ti0
         vfmsub231pd  tr1, er1, tr0
         vfmadd231pd  ti1, er1, ti0
         vmovapd      ti0, ei1
         vmovapd      tr0, er1
         vmulpd       ti2, ei2, tr0
         vmulpd       tr2, ei2, ti0
         vfmsub231pd  tr2, er2, tr0
         vfmadd231pd  ti2, er2, ti0
         vmovapd      ti0, ei2
         vmovapd      tr0, er2
         vmulpd       ti3, ei3, tr0
         vmulpd       tr3, ei3, ti0
         vfmsub231pd  tr3, er3, tr0
         vfmadd231pd  ti3, er3, ti0
         vmovapd      ti0, ei3
         vmovapd      tr0, er3
         vaddpd       er2, er0, tr0
         vaddpd       ei2, ei0, ti0
         vsubpd       er2, er0, tr2
         vsubpd       ei2, ei0, ti2
         vaddpd       er3, er1, tr1
         vaddpd       ei3, ei1, ti1
         vsubpd       er1, er3, tr3
         vsubpd       ei3, ei1, ti3
         vaddpd       tr1, tr0, er0
         vaddpd       ti1, ti0, ei0
         vaddpd       ti3, tr2, er1
         vaddpd       tr3, ti2, ei1
         vsubpd       tr1, tr0, er2
         vsubpd       ti0, ti1, ei2
         vsubpd       ti3, tr2, er3
         vsubpd       ti2, tr3, ei3
         vmovapd      er0, (%rax)
         vmovapd      er1, (%rbx)
         vmovapd      ei2, (%rcx)
         vmovapd      ei3, (%rdx)
         vmovapd      er3, (%rax, %rdi, 8)
         vmovapd      er2, (%rbx, %rdi, 8)
         vmovapd      ei1, (%rcx, %rdi, 8)
         vmovapd      ei0, (%rdx, %rdi, 8)
khiloop_cont:
         inc          khi
         cmp          N2oN1o2, khi
         jne          khiloop
         inc          klo
         cmp          N1, klo
         jne          kloloop
done:
         inc          ihi
         cmp          ihi, NHIoN1
         jne          ihiloop
         pop          %r14
         pop          %r13
         pop          %r12
         pop          %rbx
         leave
         ret
