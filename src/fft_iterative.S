#define M_SQRT         0.70710678118
#define N1 	           $4
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
	   bsf             N, %rcx
	   andq            $1, %rcx
	   jz              skipradix2
       movq            $2, %rax
       movq            N2, %rbx
       mul             %rbx
       movq            %rax, %rbx
       movq            N, %rax
       div             %rbx
       movq            %rax, NHI
       xorq            nhi, nhi
radix2loop0:
       xorq            mm, mm
radix2loop1:
	   movq			   MN2, %rax
	   shl             %rax
	   mul			   nhi
	   add			   mm, %rax
	   movq            %rax, %rcx
       lea             (X, %rax, 8), %rax
       push            %rax
       vmovapd         (%rax), er0
       lea             (%rax, MN2, 8), %rax
       push            %rax
       vmovapd         (%rax), er1
	   vaddpd          er1, er0, tr0
	   vsubpd          er1, er0, tr1
       pop             %rax
       vmovapd         tr1, (%rax)
       pop             %rbx
       vmovapd         tr0, (%rbx)
       add             SIMD_SIZE, mm
       cmp             M, mm
       jne             radix2loop1
       inc             nhi
       cmp             NHI, nhi
       jne             radix2loop0
       shl             N2
       shl             MN2
       shr             $2, NHI
       jmp             mainloop
skipradix2:
       movq            N1, %rax
       movq            N2, %rbx
       mul             %rbx
       movq            %rax, %rbx
       movq            N, %rax
       div             %rbx
       movq            %rax, NHI
mainloop:
       cmp             N, N2
       je              done
       xorq            nhi, nhi
nhiloop:
       xorq            mm, mm
mmloop0:      
	   movq			   MN2, %rax
	   shl             $2, %rax
	   mul			   nhi
	   add			   mm, %rax
	   movq            %rax, %rcx
       lea             (X, %rax, 8), %rax
       push            %rax
       vmovapd         (%rax), er0
       lea             (%rax, MN2, 8), %rax
       push            %rax
       vmovapd         (%rax), er2
       lea             (%rax, MN2, 8), %rax
       push            %rax
       vmovapd         (%rax), er1
       lea             (%rax, MN2, 8), %rax
       push            %rax
       vmovapd         (%rax), er3
	   vaddpd          er2, er0, tr0
	   vsubpd          er2, er0, tr1
	   vaddpd          er3, er1, tr2
	   vsubpd          er1, er3, tr3
	   vaddpd          tr2, tr0, er0
	   vsubpd          tr2, tr0, er1
       pop             %rax
       vmovapd         tr3, (%rax)
       pop             %rbx
       vmovapd         er1, (%rbx)
       pop             %rax
       vmovapd         tr1, (%rax)
       pop             %rbx
       vmovapd         er0, (%rbx)
       cmp             $2, N2
       jl              skipm1
	   movq            N2, %rax
	   shr             %rax
	   mul             M
	   add             %rcx, %rax
       lea             (X, %rax, 8), %rax
       push            %rax
       vmovapd         (%rax), ti0
       lea             (%rax, MN2, 8), %rax
       push            %rax
       vmovapd         (%rax), ti3
       lea             (%rax, MN2, 8), %rax
       push            %rax
       vmovapd         (%rax), ti1
       lea             (%rax, MN2, 8), %rax
       push            %rax
       vmovapd         (%rax), ti2
	   vsubpd          ti2, ti1, er1
	   vaddpd          ti2, ti1, er2
	   vmulpd          C0, er1, ti1
	   vmulpd          C1, er2, ti2
	   vaddpd          ti1, ti0, ei0
	   vsubpd          ti3, ti2, ei1
	   vsubpd          ti1, ti0, ei2
	   vaddpd          ti2, ti3, ei3
       pop             %rax
       vmovapd         ei1, (%rax)
       pop             %rax
       vmovapd         ei3, (%rax)
       pop             %rax
       vmovapd         ei2, (%rax)
       pop             %rax
       vmovapd         ei0, (%rax)
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
       mov             $3, %rcx
       mul             %rcx
       movq            %rax, %rcx
       movq            NHI, %rax
       mul             k2
       mov             $2, %rbx
       mul             %rbx
       movq            %rax, %rbx
       movq            NHI, %rax
       mul             k2
       shl             %rax
       shl             %rbx
       shl             %rcx
       vbroadcastsd    (W, %rax, 8), tr1
       vbroadcastsd    8(W, %rax, 8), ti1
       vbroadcastsd    (W, %rbx, 8), tr2
       vbroadcastsd    8(W, %rbx, 8), ti2
       vbroadcastsd    (W, %rcx, 8), tr3
       vbroadcastsd    8(W, %rcx, 8), ti3
       vmovapd         tr1, cos1
       vmovapd         ti1, sin1
       vmovapd         tr2, cos2
       vmovapd         ti2, sin2
       vmovapd         tr3, cos3
       vmovapd         ti3, sin3
       xorq            mm, mm
mmloop2:
	   movq			   MN2, %rax
	   shl             $2, %rax
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
       vmovapd         (%rax), er2
       lea             (%rax, MN2, 8), %rax
       push            %rax
       vmovapd         (%rax), er1
       lea             (%rax, MN2, 8), %rax
       push            %rax
       vmovapd         (%rax), er3
	   movq            N2, %rax
	   sub             k2, %rax
	   mul             M
	   add             %rcx, %rax
	   lea             (X, %rax, 8), %rax
       push            %rax
       vmovapd         (%rax), ei0
       lea             (%rax, MN2, 8), %rax
       push            %rax
       vmovapd         (%rax), ei2
       lea             (%rax, MN2, 8), %rax
       push            %rax
       vmovapd         (%rax), ei1
       lea             (%rax, MN2, 8), %rax
       push            %rax
       vmovapd         (%rax), ei3
       vmulpd          sin1, ei1, tr1
       vmulpd          cos1, ei1, ti1
       vfmsub231pd     cos1, er1, tr1
       vfmadd231pd     sin1, er1, ti1
       vmovapd         ti1, ei1
       vmovapd         tr1, er1
       vmulpd          sin2, ei2, tr2
       vmulpd          cos2, ei2, ti2
       vfmsub231pd     cos2, er2, tr2
       vfmadd231pd     sin2, er2, ti2
       vmovapd         ti2, ei2
       vmovapd         tr2, er2
       vmulpd          sin3, ei3, tr3
       vmulpd          cos3, ei3, ti3
       vfmsub231pd     cos3, er3, tr3
       vfmadd231pd     sin3, er3, ti3
       vmovapd         ti3, ei3
       vmovapd         tr3, er3
	   vaddpd          er2, er0, tr0
	   vaddpd          ei2, ei0, ti0
	   vsubpd          er2, er0, tr2
	   vsubpd          ei2, ei0, ti2
	   vaddpd          er3, er1, tr1
	   vaddpd          ei3, ei1, ti1
	   vsubpd          er1, er3, tr3
	   vsubpd          ei3, ei1, ti3
	   vaddpd          tr1, tr0, er0
	   vaddpd          ti1, ti0, ei0
	   vaddpd          ti3, tr2, er1
	   vaddpd          tr3, ti2, ei1
	   vsubpd          tr1, tr0, er2
	   vsubpd          ti0, ti1, ei2
	   vsubpd          ti3, tr2, er3
	   vsubpd          ti2, tr3, ei3
	   pop             %rax
       vmovapd         ei0, (%rax)
       pop             %rax
       vmovapd         ei1, (%rax)
       pop             %rax
       vmovapd         er2, (%rax)
       pop             %rax
       vmovapd         er3, (%rax)
       pop             %rax
       vmovapd         ei3, (%rax)
       pop             %rax
       vmovapd         ei2, (%rax)
       pop             %rax
       vmovapd         er1, (%rax)
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
       shl             $2, N2
       shl             $2, MN2
       shr             $2, NHI
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

       .data
none:  .double        -1.00
       .double        -1.00
       .double        -1.00
       .double        -1.00
       .double        -1.00
       .double        -1.00
       .double        -1.00
       .double        -1.00
zero:  .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
cos1:  .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
sin1:  .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
cos2:  .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
sin2:  .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
cos3:  .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
sin3:  .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
C0:    .double         M_SQRT
       .double         M_SQRT
       .double         M_SQRT
       .double         M_SQRT
       .double         M_SQRT
       .double         M_SQRT
       .double         M_SQRT
       .double         M_SQRT
C1:    .double         -M_SQRT
       .double         -M_SQRT
       .double         -M_SQRT
       .double         -M_SQRT
       .double         -M_SQRT
       .double         -M_SQRT
       .double         -M_SQRT
       .double         -M_SQRT


