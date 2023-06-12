#define  STACK_SIZE   $64
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
#define  er1          %ymm1
#define  er2          %ymm2
#define  er3          %ymm3
#define  tr0          %ymm4
#define  tr1          %ymm5
#define  tr2          %ymm6
#define  tr3          %ymm7
#define  ei0          %ymm8
#define  ei1          %ymm9
#define  ei2          %ymm10
#define  ei3          %ymm11
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
#define  khi0         khi
#define  klo0         klo

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
         mov          %rcx, N2oN1
         xor          ihi, ihi
ihiloop0:
         mov          ihi, %rax
         mov          N2N1, %rbx
         mul          %rbx
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
         inc          ihi
         cmp          ihi, NHIoN1
         jne          ihiloop0
         cmp          $1, N2
         jle          done
         mov          N2o2, %rax
         mov          %rax, %rbx
         shr          $2, %rax
         and          $3, %rbx
         mov          %rax, khi0
         mov          %rbx, klo0
         xor          ihi, ihi
ihiloop1:
         mov          NHIoN1, %rcx
         mov          N2, %rbx
         mov          klo0, %rax
         mul          %rcx
         add          ihi, %rax
         mul          %rbx
         add          khi0, %rax
         shl          $2, %rax
         lea          (X, %rax, 8), %rax
         lea          (%rax, N2, 8), %rbx
         lea          (%rbx, N2, 8), %rcx
         lea          (%rcx, N2, 8), %rdx
         vmovupd      (%rax), er0
         vmovupd      (%rbx), er1
         vmovupd      (%rcx), er2
         vmovupd      (%rdx), er3
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
         inc          ihi
         cmp          ihi, NHIoN1
         jne          ihiloop1
         cmp          $2, N2
         jle          done
         xor          ihi, ihi
ihiloop2:
         xor          klo, klo
kloloop:
         xor          khi, khi
khiloop:
         mov          khi, k2
         shl          $2, k2
         or           klo, k2
         cmp          k2, N2o2
         jle          khiloop_cont
         cmp          $0, k2
         je           khiloop_cont
         mov          N2, %rax
         sub          k2, %rax
         mov          %rax, %rbx
         shr          $2, %rax
         and          $3, %rbx
         mov          %rax, khi1
         mov          %rbx, klo1
         mov          NHI, %rax
         mov          k2, %rbx
         mul          %rbx
         shl          %rax
         mov          %rax, %rbx
         mov          %rax, %rcx
         add          %rax, %rbx
         add          %rbx, %rcx
         vbroadcastsd 0(W, %rax, 8), tr1
         vbroadcastsd 8(W, %rax, 8), ti1
         vbroadcastsd 0(W, %rbx, 8), tr2
         vbroadcastsd 8(W, %rbx, 8), ti2
         vbroadcastsd 0(W, %rcx, 8), tr3
         vbroadcastsd 8(W, %rcx, 8), ti3
         mov          NHIoN1, %rcx
         mov          N2, %rbx
         mov          klo1, %rax
         mul          %rcx
         add          ihi, %rax
         mul          %rbx
         add          khi1, %rax
         shl          $2, %rax
         mov          %rax, %rdi
         mov          klo0, %rax
         mul          %rcx
         add          ihi, %rax
         mul          %rbx
         add          khi0, %rax
         shl          $2, %rax
         sub          %rax, %rdi
         lea          (X, %rax, 8), %rax
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
         cmp          N2oN1, khi
         jne          khiloop
         inc          klo
         cmp          N1, klo
         jne          kloloop
done:
         inc          ihi
         cmp          ihi, NHIoN1
         jne          ihiloop2
         pop          %r14
         pop          %r13
         pop          %r12
         pop          %rbx
         leave
         ret
