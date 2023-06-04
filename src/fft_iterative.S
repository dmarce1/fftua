#define M_SQRT2_1      0.70710678118
#define N1 	           $4
#define SIMD_SIZE      $4
#define mm             %r10
#define nhi            %r11
#define k2             %r12
#define MN2            %r13
#define W              %r14
#define X              %r15
#define N              -32(%rbp)
#define cos1           -64(%rbp)
#define cos2           -96(%rbp)
#define cos3           -128(%rbp)
#define sin1           -160(%rbp)
#define sin2           -192(%rbp)
#define sin3           -224(%rbp)
#define NHI            -232(%rbp)
#define N2             -240(%rbp)
#define MN21           -248(%rbp)
#define MN20           -256(%rbp)
#define tmp0           -264(%rbp)
#define tmp1           -272(%rbp)
#define tmp2           -280(%rbp)
#define tmp3           -288(%rbp)
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
	   push            %rbp
	   sub             $16, %rsp
	   mov             %rsp, %rbp
	   sub             $320, %rsp
       push            %rbx
       push            %r12
       push            %r13
       push            %r14
       push            %r15
       movq            %rdi, X
       movq            %rsi, W
       movq            %rdx, N
       movq            %rcx, N
       movq            $1, N2
       xorq            %rdx, %rdx
       movq            N, %rax
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
       cmp             N, mm
       jne             radix2loop1
       inc             nhi
       cmp             NHI, nhi
       jne             radix2loop0
       mov             N2, %rax
       mov             NHI, %rbx
       shl             %rax
       shr             $2, %rbx
       mov             %rax, N2
       mov             %rbx, NHI
       shl             MN2
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
       mov             N2, %rax
       cmp             N, %rax
       je              done
       xorq            nhi, nhi
nhiloop:
       xorq            mm, mm
mmloop0:      
	   movq			   MN2, %rax
	   mul			   nhi
	   lea             (mm, %rax, 4), %rax
       lea             (X, %rax, 8), %r8
       lea             (%r8, MN2, 8), %r9
       lea             (%r9, MN2, 8), %rcx
       lea             (%rcx, MN2, 8), %rdi
       vmovapd         (%r8), er0
       vmovapd         (%r9), er2
       vmovapd         (%rcx), er1
       vmovapd         (%rdi), er3
	   vaddpd          er2, er0, tr0
	   vsubpd          er2, er0, tr1
	   vaddpd          er3, er1, tr2
	   vsubpd          er1, er3, tr3
	   vaddpd          tr2, tr0, er0
	   vsubpd          tr2, tr0, er1
       vmovapd         tr3, (%rdi)
       vmovapd         er1, (%rcx)
       vmovapd         tr1, (%r9)
       vmovapd         er0, (%r8)
       cmp             $2, N2
       jl              skipm1
	   movq            N2, %rax
	   shr             %rax
	   mulq            N
       lea             (%r8, %rax, 8), %r8
       lea             (%r9, %rax, 8), %r9
       lea             (%rcx, %rax, 8), %rcx
       lea             (%rdi, %rax, 8), %rdi
       vmovapd         (%r8), ti0
       vmovapd         (%r9), ti3
       vmovapd         (%rcx), ti1
       vmovapd         (%rdi), ti2
	   vsubpd          ti2, ti1, er1
	   vaddpd          ti2, ti1, er2
	   vmulpd          C0, er1, ti1
	   vmulpd          C1, er2, ti2
	   vaddpd          ti1, ti0, ei0
	   vsubpd          ti3, ti2, ei1
	   vsubpd          ti1, ti0, ei2
	   vaddpd          ti2, ti3, ei3
       vmovapd         ei1, (%rdi)
       vmovapd         ei3, (%rcx)
       vmovapd         ei2, (%r9)
       vmovapd         ei0, (%r8)
skipm1:
       add             SIMD_SIZE, mm
       cmp             N, mm
       jne             mmloop0
       cmp             $2, N2
       jle             skipk2
       movq            $1, k2
k2loop:
       movq            NHI, %rax
       mul             k2
       shl             %rax
       mov             %rax, %rcx
       add             %rcx, %rax
       mov             %rax, %rbx
       add             %rcx, %rax
       vbroadcastsd    (W, %rcx, 8), tr1
       vbroadcastsd    8(W, %rcx, 8), ti1
       vbroadcastsd    (W, %rbx, 8), tr2
       vbroadcastsd    8(W, %rbx, 8), ti2
       vbroadcastsd    (W, %rax, 8), tr3
       vbroadcastsd    8(W, %rax, 8), ti3
       vmovapd         tr1, cos1
       vmovapd         ti1, sin1
       vmovapd         tr2, cos2
       vmovapd         ti2, sin2
       vmovapd         tr3, cos3
       vmovapd         ti3, sin3
       xorq            mm, mm
mmloop2:
 	   movq			   MN2, %rax
	   mul			   nhi
	   lea             (mm, %rax, 4), %rcx
	   movq            N2, %rax
	   sub             k2, %rax
	   sub             k2, %rax
	   mulq            N
	   mov             %rax, %rbx
	   movq            N, %rax
	   mul             k2
	   add             %rcx, %rax
       lea             (X, %rax, 8), %r8
       lea             (%r8, MN2, 8), %r9
       lea             (%r9, MN2, 8), %rcx
       lea             (%rcx, MN2, 8), %rdx
       vmovapd         (%r8), er0
       vmovapd         (%r9), er2
       vmovapd         (%rcx), er1
       vmovapd         (%rdx), er3
       vmovapd         (%r8, %rbx, 8), ei0
       vmovapd         (%r9, %rbx, 8), ei2
       vmovapd         (%rcx, %rbx, 8), ei1
       vmovapd         (%rdx, %rbx, 8), ei3
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
       vmovapd         ei3, (%rdx)
       vmovapd         ei2, (%rcx)
       vmovapd         er1, (%r9)
       vmovapd         er0, (%r8)
       vmovapd         ei0, (%rdx, %rbx, 8)
       vmovapd         ei1, (%rcx, %rbx, 8)
       vmovapd         er2, (%r9, %rbx, 8)
       vmovapd         er3, (%r8, %rbx, 8)
       add             SIMD_SIZE, mm
       cmp             N, mm
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
       mov             N2, %rax
       mov             NHI, %rbx
       shl             $2, %rax
       shr             $2, %rbx
       mov             %rax, N2
       mov             %rbx, NHI
       shl             $2, MN2
       jmp             mainloop
done:
       pop             %r15
       pop             %r14
       pop             %r13
       pop             %r12
       pop             %rbx
       emms
 	   add             $336, %rsp
       pop             %rbp
       ret

       .data
none:  .double        -1.00
       .double        -1.00
       .double        -1.00
       .double        -1.00
zero:  .double         0.00
       .double         0.00
       .double         0.00
       .double         0.00
C0:    .double         M_SQRT2_1
       .double         M_SQRT2_1
       .double         M_SQRT2_1
       .double         M_SQRT2_1
C1:    .double         -M_SQRT2_1
       .double         -M_SQRT2_1
       .double         -M_SQRT2_1
       .double         -M_SQRT2_1
two:   .quad           2
       .quad           2
       .quad           2
       .quad           2

