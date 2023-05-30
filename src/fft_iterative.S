#define N1 	      $2
#define SIMD_SIZE $4
#define mm        %r8
#define M         %r9
#define nhi       %r10
#define NHI       %r11
#define k2        %r12
#define N2        %r13
#define MN2       %r14
#define W         %r15
#define X         %rdi
#define N         %rsi
#define MN2N1     %rcx
#define er0       %ymm0
#define er1       %ymm1
#define ei0       %ymm2
#define ei1       %ymm3
#define tr0       %ymm4
#define tr1       %ymm5
#define ti0       %ymm6
#define ti1       %ymm7
#define two       %ymm8
#define zero      %ymm9

        .global        fft_iter_real

       .data
dtwo:  .double         2.00
dzero: .double         0.00

       .text
fft_iter_real:
       push            %rbx
       push            %rdi
       push            %rsi
       push            %r8
       push            %r9
       push            %r10
       push            %r11
       push            %r12
       push            %r13
       push            %r14
       push            %r15
       movq            %rdi, X
       movq            %rsi, W
       movq            %rdx, N
       movq            %rcx, M
       movq            $1, N2
       vbroadcastsd    dtwo, two
       vbroadcastsd    dzero, zero
mainloop:
       movq            N2, %rax
       cmp             N, %rax
       je              done
       xorq            %rdx, %rdx
       mov             M, %rax
       mov             N2, MN2
       mul             MN2
       movq            %rax, MN2
       movq            N1, %rax
       movq            N2, %rbx
       mul             %rbx
       movq            %rax, %rbx
       movq            N, %rax
       div             %rbx
       movq            %rax, NHI
       movq			   N1, %rax
       movq			   MN2, MN2N1
       mul             MN2N1
       mov		       %rax, MN2N1
       xorq            nhi, nhi
nhiloop:
       xorq            mm, mm
mmloop0:      
	   movq            MN2N1, %rax
	   mul             nhi
       lea             (X, %rax, 8), %rax
       lea             (%rax, mm, 8), %rax
       push            %rax
       vmovapd         (%rax), er0
       lea             (%rax, MN2, 8), %rax
       push            %rax
       vmovapd         (%rax), er1
       vmovapd         er1, tr1
       vsubpd          er1, er0, er1
       vaddpd          tr1, er0, er0
       pop             %rax
       vmovapd         er1, (%rax)
       pop             %rax
       vmovapd         er0, (%rax)
       add             SIMD_SIZE, mm
       cmp             M, mm
       jne             mmloop0
       cmp             $2, N2
       jl              skipk2
       xorq            mm, mm
mmloop1:
       movq            MN2, %rax
       shr             %rax
       add             MN2, %rax
       push            %rax
	   movq            MN2N1, %rax
	   mul             nhi
	   pop             %rbx
	   add			   %rbx, %rax
       add             mm, %rax
       lea             (X, %rax, 8), %rax
       vmovapd         (%rax), er1
       vsubpd          er1, zero, er1
       vmovapd         er1, (%rax)
       add             SIMD_SIZE, mm
       cmp             M, mm
       jne             mmloop1
       movq            $1, k2
k2loop:      
       cmp             $2, N2
       jle             skipk2
       movq            NHI, %rax
       xorq            %rdx, %rdx
       mul             k2
       vbroadcastsd    (W, %rax, 8), tr1
       vbroadcastsd    8(W, %rax, 8), ti1
       xorq            mm, mm
mmloop2:
       mov             M, %rax
       xorq            %rdx, %rdx
       mul             k2
       lea             (X, %rax, 8), %rax
       lea             (%rax, mm, 8), %rax
       push            %rax
       vmovapd         (%rax), er0
       add             MN2, %rax
       push            %rax
       vmovapd         (%rax), er1
       mov             N2, %rax
       sub             k2, %rax
       xorq            %rdx, %rdx
       mul             M
       lea             (X, %rax, 8), %rax
       lea             (%rax, mm, 8), %rax
       push            %rax
       vmovapd         (%rax), ei0
       lea             (%rax, MN2, 8), %rax
       push            %rax
       vmovapd         (%rax), ei1
       vmovapd         er1, tr0
       vmovapd         ei1, ti0
       vfmsub132pd     ti0, er0, ti1
       vfmsub132pd     tr0, er0, tr1
       vfmadd132pd     ti0, ei0, tr1
       vfmadd132pd     tr0, ei0, ti1
       vmovapd         er0, er1
       vmovapd         ei0, ei1
       vfmsub132pd     tr0, er1, two
       vfmsub132pd     ti0, ei1, two
       vsubpd          ei1, zero, ei1
       pop             %rax
       vmovapd         ei0, (%rax)
       pop             %rax
       vmovapd         ei1, (%rax)
       pop             %rax
       vmovapd         er1, (%rax)
       pop             %rax
       vmovapd         er0, (%rax)
       add             SIMD_SIZE, mm
       cmp             M, mm
       jne             mmloop2
       inc             k2
       mov             N2, %rax
       shr             %rax
       cmp             k2, %rax
       jne             k2loop
skipk2:
       inc             nhi
       cmp             NHI, nhi
       jne             nhiloop
       shl             N2
       jmp             mainloop
done:
       pop             %r15
       pop             %r14
       pop             %r13
       pop             %r12
       pop             %r11
       pop             %r10
       pop             %r9
       pop             %r8
       pop             %rsi
       pop             %rdi
       pop             %rbx
       ret


