#define N1 	           $2
#define SIMD_SIZE      $4
#define mm             %r8
#define M              %r9
#define nhi            %r10
#define NHI            %r11
#define k2             %r12
#define N2             %r13
#define MN2            %r14
#define W              %r15
#define X              %rdi
#define N              %rsi
#define er0            %ymm0
#define er1            %ymm1
#define er2            %ymm2
#define er3            %ymm3
#define ei0            %ymm4
#define ei1            %ymm5
#define ei2            %ymm6
#define ei3            %ymm7
#define tr0            %ymm8
#define tr1            %ymm9
#define tr2            %ymm10
#define tr3            %ymm11
#define ti0            %ymm12
#define ti1            %ymm13
#define ti2            %ymm14
#define ti3            %ymm15

       .global         fft_iter_real

       .data
two:   .double         2.00
       .double         2.00
       .double         2.00
       .double         2.00
       .double         2.00
       .double         2.00
       .double         2.00
       .double         2.00
zero:  .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
wr:    .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
wi:    .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00

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
       xorq            %rdx, %rdx
       movq            M, %rax
       movq            N2, MN2
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
mainloop:
       cmp             N, N2
       je              done
       xorq            nhi, nhi
nhiloop:
       xorq            mm, mm
mmloop0:      
	   movq			   MN2, %rax
	   shl             %rax
	   mul			   nhi
	   add			   mm, %rax
	   mov             %rax, %rcx
	   lea             (X, %rax, 8), %rax
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
       cmp             $2, N2
       jl              skipm1
	   movq            N2, %rax
	   shr             %rax
	   add			   N2, %rax
	   mul             M
	   add             %rcx, %rax
       lea             (X, %rax, 8), %rax
       vmovapd         zero, tr1
       vsubpd          (%rax), tr1, er1
       vmovapd         er1, (%rax)
skipm1:
       add             SIMD_SIZE, mm
       cmp             M, mm
       jne             mmloop0
       cmp             $2, N2
       jle             skipk2
       movq            $1, k2
k2loop:
       movq            NHI, %rax
       mul             k2
       shl             %rax
       vbroadcastsd    (W, %rax, 8), tr0
       vbroadcastsd    8(W, %rax, 8), ti0
       vmovapd         tr0, wr
       vmovapd         ti0, wi
       xorq            mm, mm
mmloop2:
	   movq			   MN2, %rax
	   shl             %rax
	   mul			   nhi
	   add			   mm, %rax
	   mov             %rax, %rcx
	   movq            M, %rax
	   mul             k2
	   add             %rcx, %rax
       lea             (X, %rax, 8), %rax
       push            %rax
       vmovapd         (%rax), er0
       lea             (%rax, MN2, 8), %rax
       push            %rax
       vmovapd         (%rax), er1
	   movq            N2, %rax
	   sub             k2, %rax
	   mul             M
	   add             %rcx, %rax
	   lea             (X, %rax, 8), %rax
       push            %rax
       vmovapd         (%rax), ei0
       lea             (%rax, MN2, 8), %rax
       push            %rax
       vmovapd         (%rax), ei1
       vmovapd         er0, tr0
       vmovapd         ei0, ti0
       vfmsub231pd     wi, ei1, er0
       vfmsub231pd     wr, er1, er0
       vfmadd231pd     wr, ei1, ei0
       vfmadd231pd     wi, er1, ei0
       vmovapd		   er0, er1
       vmovapd		   ei0, ei1
       vfmsub231pd     two, tr0, er1
       vfnmadd231pd    two, ti0, ei1
       pop             %rax
       vmovapd         ei0, (%rax)
       pop             %rax
       vmovapd         er1, (%rax)
       pop             %rax
       vmovapd         ei1, (%rax)
       pop             %rax
       vmovapd         er0, (%rax)
       add             SIMD_SIZE, mm
       cmp             M, mm
       jne             mmloop2
       inc             k2
       movq            N2, %rax
       shr             %rax
       cmp             k2, %rax
       jne             k2loop
skipk2:
       inc             nhi
       cmp             NHI, nhi
       jne             nhiloop
       shl             N2
       shl             MN2
       shr             NHI
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
       emms
       ret


