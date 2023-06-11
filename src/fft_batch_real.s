	.file	"fft_batch_real.cpp"
# GNU C++14 (Ubuntu 9.4.0-1ubuntu1~20.04.1) version 9.4.0 (x86_64-linux-gnu)
#	compiled by GNU C version 9.4.0, GMP version 6.2.0, MPFR version 4.0.2, MPC version 1.1.0, isl version isl-0.22.1-GMP

# GGC heuristics: --param ggc-min-expand=100 --param ggc-min-heapsize=131072
# options passed:  -I ../include/ -imultiarch x86_64-linux-gnu
# -D_GNU_SOURCE fft_batch_real.cpp -march=znver2 -mmmx -mno-3dnow -msse
# -msse2 -msse3 -mssse3 -msse4a -mcx16 -msahf -mmovbe -maes -msha -mpclmul
# -mpopcnt -mabm -mno-lwp -mfma -mno-fma4 -mno-xop -mbmi -mno-sgx -mbmi2
# -mno-pconfig -mwbnoinvd -mno-tbm -mavx -mavx2 -msse4.2 -msse4.1 -mlzcnt
# -mno-rtm -mno-hle -mrdrnd -mf16c -mfsgsbase -mrdseed -mprfchw -madx
# -mfxsr -mxsave -mxsaveopt -mno-avx512f -mno-avx512er -mno-avx512cd
# -mno-avx512pf -mno-prefetchwt1 -mclflushopt -mxsavec -mxsaves
# -mno-avx512dq -mno-avx512bw -mno-avx512vl -mno-avx512ifma -mno-avx512vbmi
# -mno-avx5124fmaps -mno-avx5124vnniw -mclwb -mmwaitx -mclzero -mno-pku
# -mrdpid -mno-gfni -mno-shstk -mno-avx512vbmi2 -mno-avx512vnni -mno-vaes
# -mno-vpclmulqdq -mno-avx512bitalg -mno-avx512vpopcntdq -mno-movdiri
# -mno-movdir64b -mno-waitpkg -mno-cldemote -mno-ptwrite
# --param l1-cache-size=32 --param l1-cache-line-size=64
# --param l2-cache-size=512 -mtune=znver2 -O3 -fverbose-asm
# -fasynchronous-unwind-tables -fstack-protector-strong -Wformat
# -Wformat-security -fstack-clash-protection -fcf-protection
# options enabled:  -fPIC -fPIE -faggressive-loop-optimizations
# -falign-functions -falign-jumps -falign-labels -falign-loops
# -fassume-phsa -fasynchronous-unwind-tables -fauto-inc-dec
# -fbranch-count-reg -fcaller-saves -fcode-hoisting
# -fcombine-stack-adjustments -fcommon -fcompare-elim -fcprop-registers
# -fcrossjumping -fcse-follow-jumps -fdefer-pop
# -fdelete-null-pointer-checks -fdevirtualize -fdevirtualize-speculatively
# -fdwarf2-cfi-asm -fearly-inlining -feliminate-unused-debug-types
# -fexceptions -fexpensive-optimizations -fforward-propagate
# -ffp-int-builtin-inexact -ffunction-cse -fgcse -fgcse-after-reload
# -fgcse-lm -fgnu-runtime -fgnu-unique -fguess-branch-probability
# -fhoist-adjacent-loads -fident -fif-conversion -fif-conversion2
# -findirect-inlining -finline -finline-atomics -finline-functions
# -finline-functions-called-once -finline-small-functions -fipa-bit-cp
# -fipa-cp -fipa-cp-clone -fipa-icf -fipa-icf-functions -fipa-icf-variables
# -fipa-profile -fipa-pure-const -fipa-ra -fipa-reference
# -fipa-reference-addressable -fipa-sra -fipa-stack-alignment -fipa-vrp
# -fira-hoist-pressure -fira-share-save-slots -fira-share-spill-slots
# -fisolate-erroneous-paths-dereference -fivopts -fkeep-static-consts
# -fleading-underscore -flifetime-dse -floop-interchange
# -floop-unroll-and-jam -flra-remat -flto-odr-type-merging -fmath-errno
# -fmerge-constants -fmerge-debug-strings -fmove-loop-invariants
# -fomit-frame-pointer -foptimize-sibling-calls -foptimize-strlen
# -fpartial-inlining -fpeel-loops -fpeephole -fpeephole2 -fplt
# -fpredictive-commoning -fprefetch-loop-arrays -free -freg-struct-return
# -freorder-blocks -freorder-blocks-and-partition -freorder-functions
# -frerun-cse-after-loop -fsched-critical-path-heuristic
# -fsched-dep-count-heuristic -fsched-group-heuristic -fsched-interblock
# -fsched-last-insn-heuristic -fsched-rank-heuristic -fsched-spec
# -fsched-spec-insn-heuristic -fsched-stalled-insns-dep -fschedule-fusion
# -fschedule-insns2 -fsemantic-interposition -fshow-column -fshrink-wrap
# -fshrink-wrap-separate -fsigned-zeros -fsplit-ivs-in-unroller
# -fsplit-loops -fsplit-paths -fsplit-wide-types -fssa-backprop
# -fssa-phiopt -fstack-clash-protection -fstack-protector-strong
# -fstdarg-opt -fstore-merging -fstrict-aliasing
# -fstrict-volatile-bitfields -fsync-libcalls -fthread-jumps
# -ftoplevel-reorder -ftrapping-math -ftree-bit-ccp -ftree-builtin-call-dce
# -ftree-ccp -ftree-ch -ftree-coalesce-vars -ftree-copy-prop -ftree-cselim
# -ftree-dce -ftree-dominator-opts -ftree-dse -ftree-forwprop -ftree-fre
# -ftree-loop-distribute-patterns -ftree-loop-distribution
# -ftree-loop-if-convert -ftree-loop-im -ftree-loop-ivcanon
# -ftree-loop-optimize -ftree-loop-vectorize -ftree-parallelize-loops=
# -ftree-partial-pre -ftree-phiprop -ftree-pre -ftree-pta -ftree-reassoc
# -ftree-scev-cprop -ftree-sink -ftree-slp-vectorize -ftree-slsr -ftree-sra
# -ftree-switch-conversion -ftree-tail-merge -ftree-ter -ftree-vrp
# -funit-at-a-time -funswitch-loops -funwind-tables -fverbose-asm
# -fversion-loops-for-strides -fzero-initialized-in-bss
# -m128bit-long-double -m64 -m80387 -mabm -madx -maes -malign-stringops
# -mavx -mavx2 -mbmi -mbmi2 -mclflushopt -mclwb -mclzero -mcx16 -mf16c
# -mfancy-math-387 -mfma -mfp-ret-in-387 -mfsgsbase -mfxsr -mglibc
# -mieee-fp -mlong-double-80 -mlzcnt -mmmx -mmovbe -mmwaitx -mpclmul
# -mpopcnt -mprfchw -mpush-args -mrdpid -mrdrnd -mrdseed -mred-zone -msahf
# -msha -msse -msse2 -msse3 -msse4 -msse4.1 -msse4.2 -msse4a -mssse3 -mstv
# -mtls-direct-seg-refs -mvzeroupper -mwbnoinvd -mxsave -mxsavec -mxsaveopt
# -mxsaves

	.text
	.p2align 4
	.globl	_Z14fft_batch_realPdii
	.type	_Z14fft_batch_realPdii, @function
_Z14fft_batch_realPdii:
.LFB9575:
	.cfi_startproc
	endbr64	
	leaq	8(%rsp), %r10	#,
	.cfi_def_cfa 10, 0
	andq	$-32, %rsp	#,
# fft_batch_real.cpp:5: 	int N1 = (ilogb(N) % 2 == 0) ? 4 : 8;
	vxorps	%xmm0, %xmm0, %xmm0	# tmp455
# fft_batch_real.cpp:4: void fft_batch_real(double* x, int L, int N) {
	pushq	-8(%r10)	#
	pushq	%rbp	#
# fft_batch_real.cpp:5: 	int N1 = (ilogb(N) % 2 == 0) ? 4 : 8;
	vcvtsi2sdl	%edx, %xmm0, %xmm0	# N, tmp455, tmp456
# fft_batch_real.cpp:4: void fft_batch_real(double* x, int L, int N) {
	movq	%rsp, %rbp	#,
	.cfi_escape 0x10,0x6,0x2,0x76,0
	pushq	%r15	#
	pushq	%r14	#
	.cfi_escape 0x10,0xf,0x2,0x76,0x78
	.cfi_escape 0x10,0xe,0x2,0x76,0x70
	movq	%rdi, %r15	# tmp449, x
	pushq	%r13	#
	pushq	%r12	#
	pushq	%r10	#
	.cfi_escape 0xf,0x3,0x76,0x58,0x6
	.cfi_escape 0x10,0xd,0x2,0x76,0x68
	.cfi_escape 0x10,0xc,0x2,0x76,0x60
	pushq	%rbx	#
	.cfi_escape 0x10,0x3,0x2,0x76,0x50
	movl	%edx, %ebx	# tmp451, N
	addq	$-128, %rsp	#,
# fft_batch_real.cpp:4: void fft_batch_real(double* x, int L, int N) {
	movl	%esi, -140(%rbp)	# tmp450, %sfp
	movl	%edx, -144(%rbp)	# N, %sfp
# fft_batch_real.cpp:5: 	int N1 = (ilogb(N) % 2 == 0) ? 4 : 8;
	call	ilogb@PLT	#
# fft_batch_real.cpp:5: 	int N1 = (ilogb(N) % 2 == 0) ? 4 : 8;
	testb	$1, %al	#, tmp452
	jne	.L2	#,
# fft_batch_real.cpp:7: 	int NHI = N / (N1 * N2);
	testl	%ebx, %ebx	# N
	leal	3(%rbx), %eax	#, tmp320
	movl	%ebx, %edi	# N, N
	cmovns	%ebx, %eax	# tmp320,, N, N
	sarl	$2, %eax	#, N
	movl	%eax, %r12d	# N, NHI
# fft_batch_real.cpp:8: 	const auto& W = twiddles(N);
	call	_Z8twiddlesi@PLT	#
# fft_batch_real.cpp:6: 	int N2 = 1;
	movl	$1, -60(%rbp)	#, %sfp
# fft_batch_real.cpp:8: 	const auto& W = twiddles(N);
	movq	%rax, -96(%rbp)	# tmp453, %sfp
.L3:
# fft_batch_real.cpp:79: 	while (N2 < N) {
	movl	-60(%rbp), %ebx	# %sfp, N2
	cmpl	%ebx, -144(%rbp)	# N2, %sfp
	jle	.L36	#,
	movslq	-140(%rbp), %rdi	# %sfp,
	leaq	32(%r15), %r14	#, tmp446
	movq	%rdi, %rbx	#,
	leal	-1(%rdi), %eax	#, tmp363
	movq	%rdi, -104(%rbp)	# _561, %sfp
	shrl	$2, %eax	#,
	negl	%ebx	# tmp366
	leaq	0(,%rax,4), %r13	#, _451
	movslq	%ebx, %rax	# tmp366,
	movq	%r15, %rbx	# x, tmp445
	subq	%rdi, %rax	# _561, tmp368
	salq	$3, %rax	#, _559
	movq	%rax, -112(%rbp)	# _559, %sfp
	leaq	0(,%rdi,8), %rax	#, _550
	subq	%rax, %rbx	# _550, tmp445
	movl	-60(%rbp), %eax	# %sfp, D
	movq	%rbx, -152(%rbp)	# tmp445, %sfp
.L16:
# fft_batch_real.cpp:80: 		int D = L * N2;
	movl	-140(%rbp), %edi	# %sfp, L
	imull	%edi, %eax	# L, D
	movl	%eax, -72(%rbp)	# D, %sfp
# fft_batch_real.cpp:81: 		for (int ihi = 0; ihi < NHI; ihi++) {
	testl	%r12d, %r12d	# NHI
	jle	.L10	#,
# fft_batch_real.cpp:86: 				xi1 = xi0 + D;
	movslq	%eax, %rbx	# D, _378
# fft_batch_real.cpp:86: 				xi1 = xi0 + D;
	leaq	0(,%rbx,8), %rcx	#, _17
	testl	%edi, %edi	# L
	jle	.L11	#,
	leal	0(,%rax,4), %r10d	#, tmp369
	movq	%rbx, %r8	# _378, tmp370
	xorl	%edi, %edi	# ivtmp.106
# fft_batch_real.cpp:81: 		for (int ihi = 0; ihi < NHI; ihi++) {
	xorl	%r9d, %r9d	# ihi
	movslq	%r10d, %r10	# tmp369, _441
	salq	$4, %r8	#, tmp370
	leaq	0(,%r10,8), %r11	#, _435
	addq	%r15, %r8	# x, ivtmp.108
	.p2align 4
	.p2align 3
.L14:
	leaq	(%rdi,%r13), %rdx	#, tmp372
	leaq	(%r15,%rdi,8), %rax	#, ivtmp.97
	leaq	(%r14,%rdx,8), %rsi	#, _446
	movq	%r8, %rdx	# ivtmp.108, ivtmp.98
	.p2align 4
	.p2align 3
.L12:
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:862:   return *(__m256d *)__P;
	vmovapd	(%rax), %ymm2	# MEM[base: _464, offset: 0B], _172
	vmovapd	(%rdx), %ymm3	# MEM[base: _462, offset: 0B], _173
	addq	$32, %rax	#, ivtmp.97
	addq	$32, %rdx	#, ivtmp.98
	vmovapd	-32(%rax,%rcx), %ymm5	# MEM[base: _464, index: _17, offset: 0B], _174
	vmovapd	-32(%rdx,%rcx), %ymm0	# MEM[base: _462, index: _17, offset: 0B], _175
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:143:   return (__m256d) ((__v4df)__A + (__v4df)__B);
	vaddpd	%ymm5, %ymm2, %ymm1	# _174, _172, _176
	vaddpd	%ymm0, %ymm3, %ymm4	# _175, _173, _178
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:362:   return (__m256d) ((__v4df)__A - (__v4df)__B);
	vsubpd	%ymm5, %ymm2, %ymm2	# _174, _172, tmp376
	vsubpd	%ymm3, %ymm0, %ymm0	# _173, _175, tmp378
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:143:   return (__m256d) ((__v4df)__A + (__v4df)__B);
	vaddpd	%ymm4, %ymm1, %ymm6	# _178, _176, tmp375
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:362:   return (__m256d) ((__v4df)__A - (__v4df)__B);
	vsubpd	%ymm4, %ymm1, %ymm1	# _178, _176, tmp377
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:868:   *(__m256d *)__P = __A;
	vmovapd	%ymm6, -32(%rax)	# tmp375, MEM[base: _464, offset: 0B]
	vmovapd	%ymm2, -32(%rax,%rcx)	# tmp376, MEM[base: _464, index: _17, offset: 0B]
	vmovapd	%ymm1, -32(%rdx)	# tmp377, MEM[base: _462, offset: 0B]
	vmovapd	%ymm0, -32(%rdx,%rcx)	# tmp378, MEM[base: _462, index: _17, offset: 0B]
# fft_batch_real.cpp:82: 			for (int l = 0; l < L; l += SIMD_SIZE) {
	cmpq	%rax, %rsi	# ivtmp.97, _446
	jne	.L12	#,
# fft_batch_real.cpp:81: 		for (int ihi = 0; ihi < NHI; ihi++) {
	incl	%r9d	# ihi
	addq	%r10, %rdi	# _441, ivtmp.106
	addq	%r11, %r8	# _435, ivtmp.108
# fft_batch_real.cpp:81: 		for (int ihi = 0; ihi < NHI; ihi++) {
	cmpl	%r9d, %r12d	# ihi, NHI
	jne	.L14	#,
# fft_batch_real.cpp:105: 		if (N2 > 1) {
	cmpl	$1, -60(%rbp)	#, %sfp
	jne	.L39	#,
.L10:
# fft_batch_real.cpp:205: 		N2 *= N1;
	movl	-60(%rbp), %eax	# %sfp, N2
	leal	0(,%rax,4), %ecx	#, N2
# fft_batch_real.cpp:207: 		NHI = N / (N1 * N2);
	sall	$4, %eax	#, N2
	movl	%eax, %esi	# N2, tmp379
# fft_batch_real.cpp:207: 		NHI = N / (N1 * N2);
	movl	-144(%rbp), %eax	# %sfp, NHI
	cltd
	idivl	%esi	# tmp379
	movl	%eax, %r12d	# NHI, NHI
# fft_batch_real.cpp:79: 	while (N2 < N) {
	cmpl	%ecx, -144(%rbp)	# N2, %sfp
	jle	.L36	#,
	movl	%ecx, -60(%rbp)	# N2, %sfp
	movl	%ecx, %eax	# N2, D
	jmp	.L16	#
.L39:
# fft_batch_real.cpp:108: 					int i0 = N2 / 2 + N2 * N1 * ihi;
	movl	-60(%rbp), %edi	# %sfp, tmp436
	movl	-72(%rbp), %eax	# %sfp, D
	movq	%rbx, %r10	# _378, tmp438
	xorl	%r8d, %r8d	# ihi
	salq	$4, %r10	#, tmp438
	addq	%r15, %r10	# x, _481
	sarl	%edi	# tmp436
	imull	-140(%rbp), %edi	# %sfp, tmp437
	leal	0(,%rax,4), %r9d	#, tmp435
	movslq	%r9d, %r9	# tmp435, _496
	movslq	%edi, %rdi	# tmp437, ivtmp.87
.L21:
	leaq	0(,%rdi,8), %rdx	#, _485
	leaq	0(%r13,%rdi), %rsi	#, tmp382
	leaq	(%r15,%rdx), %rax	#, ivtmp.78
	leaq	(%r14,%rsi,8), %rsi	#, _501
	addq	%r10, %rdx	# _481, ivtmp.79
	.p2align 4
	.p2align 3
.L17:
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:362:   return (__m256d) ((__v4df)__A - (__v4df)__B);
	vmovapd	(%rdx), %ymm7	# MEM[base: _523, offset: 0B], tmp632
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:143:   return (__m256d) ((__v4df)__A + (__v4df)__B);
	vmovapd	.LC1(%rip), %ymm3	#, _168
	addq	$32, %rax	#, ivtmp.78
	addq	$32, %rdx	#, ivtmp.79
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:362:   return (__m256d) ((__v4df)__A - (__v4df)__B);
	vsubpd	-32(%rdx,%rcx), %ymm7, %ymm1	# MEM[base: _523, index: _17, offset: 0B], tmp632, _161
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:143:   return (__m256d) ((__v4df)__A + (__v4df)__B);
	vaddpd	-32(%rdx,%rcx), %ymm7, %ymm0	# MEM[base: _523, index: _17, offset: 0B], tmp633, _165
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:362:   return (__m256d) ((__v4df)__A - (__v4df)__B);
	vmovapd	-32(%rax), %ymm7	# MEM[base: _525, offset: 0B], tmp636
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:143:   return (__m256d) ((__v4df)__A + (__v4df)__B);
	vmovapd	-32(%rax,%rcx), %ymm6	# MEM[base: _525, index: _17, offset: 0B], tmp637
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:362:   return (__m256d) ((__v4df)__A - (__v4df)__B);
	vmovapd	.LC0(%rip), %ymm2	#, _169
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:143:   return (__m256d) ((__v4df)__A + (__v4df)__B);
	vfmadd213pd	-32(%rax), %ymm1, %ymm3	# MEM[base: _525, offset: 0B], _161, _168
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:362:   return (__m256d) ((__v4df)__A - (__v4df)__B);
	vfnmadd132pd	.LC1(%rip), %ymm7, %ymm1	#, tmp636, _170
	vfmsub213pd	-32(%rax,%rcx), %ymm0, %ymm2	# MEM[base: _525, index: _17, offset: 0B], _165, _169
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:143:   return (__m256d) ((__v4df)__A + (__v4df)__B);
	vfmadd132pd	.LC0(%rip), %ymm6, %ymm0	#, tmp637, _171
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:868:   *(__m256d *)__P = __A;
	vmovapd	%ymm3, -32(%rax)	# _168, MEM[base: _525, offset: 0B]
	vmovapd	%ymm1, -32(%rax,%rcx)	# _170, MEM[base: _525, index: _17, offset: 0B]
	vmovapd	%ymm0, -32(%rdx)	# _171, MEM[base: _523, offset: 0B]
	vmovapd	%ymm2, -32(%rdx,%rcx)	# _169, MEM[base: _523, index: _17, offset: 0B]
# fft_batch_real.cpp:107: 				for (int l = 0; l < L; l += SIMD_SIZE) {
	cmpq	%rax, %rsi	# ivtmp.78, _501
	jne	.L17	#,
# fft_batch_real.cpp:106: 			for (int ihi = 0; ihi < NHI; ihi++) {
	incl	%r8d	# ihi
	addq	%r9, %rdi	# _496, ivtmp.87
# fft_batch_real.cpp:106: 			for (int ihi = 0; ihi < NHI; ihi++) {
	cmpl	%r8d, %r12d	# ihi, NHI
	jne	.L21	#,
.L20:
# fft_batch_real.cpp:133: 		for (int k2 = 1; k2 < N2 / 2; k2++) {
	movl	-60(%rbp), %eax	# %sfp, N2
	movl	%eax, %edi	# N2, _552
	sarl	%edi	# _552
	movl	%edi, -64(%rbp)	# _552, %sfp
# fft_batch_real.cpp:133: 		for (int k2 = 1; k2 < N2 / 2; k2++) {
	cmpl	$3, %eax	#, N2
	jle	.L10	#,
	movl	-140(%rbp), %eax	# %sfp,
# fft_batch_real.cpp:152: 					xi2 = xi0 + D;
	leaq	0(,%rbx,8), %rsi	#, _502
	testl	%eax, %eax	#
	jle	.L10	#,
	movslq	%r12d, %rax	# NHI, NHI
	movl	-72(%rbp), %edi	# %sfp, D
	movq	-152(%rbp), %r11	# %sfp, tmp445
	salq	$4, %rax	#, NHI
	movq	%rax, %rcx	# NHI, _613
	movq	%rax, -136(%rbp)	# _613, %sfp
	leal	(%r12,%r12), %eax	#, tmp390
	movslq	%eax, %rdx	# tmp390, tmp391
	addl	%r12d, %eax	# NHI, tmp394
	movq	%rcx, -88(%rbp)	# _613, %sfp
	cltq
	salq	$4, %rdx	#, tmp391
	salq	$4, %rax	#, tmp395
	movq	%rdx, -128(%rbp)	# _606, %sfp
	movq	%rdx, -80(%rbp)	# _606, %sfp
	movq	%rax, %r10	# tmp395, _595
	movq	%rax, -120(%rbp)	# _595, %sfp
	movl	-60(%rbp), %eax	# %sfp, N2
	decl	%eax	# ivtmp.69
	movl	%eax, -56(%rbp)	# ivtmp.69, %sfp
	movl	%edi, %eax	# D, tmp396
	subl	-140(%rbp), %eax	# %sfp, tmp396
	sall	$2, %edi	#, D
	cltq
	leaq	(%r11,%rax,8), %rax	#, ivtmp.72
	movslq	%edi, %r11	# D,
	leaq	(%rbx,%rbx,2), %rdi	#, tmp403
	movq	-104(%rbp), %rbx	# %sfp, _561
	salq	$3, %rdi	#, tmp404
	movq	%rbx, -72(%rbp)	# _561, %sfp
	movq	%r10, %rbx	# _595, ivtmp.66
	.p2align 4
	.p2align 3
.L25:
# fft_batch_real.cpp:135: 			const auto w1 = W[1 * j];
	movq	-96(%rbp), %rcx	# %sfp, _273
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:1306:   return __extension__ (__m256d){ __A, __A, __A, __A };
	movq	-80(%rbp), %r10	# %sfp, ivtmp.65
	movq	-72(%rbp), %r9	# %sfp, ivtmp.53
# fft_batch_real.cpp:135: 			const auto w1 = W[1 * j];
	movq	(%rcx), %rdx	# MEM[(struct complex * *)_273], _232
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:1306:   return __extension__ (__m256d){ __A, __A, __A, __A };
	movq	-88(%rbp), %rcx	# %sfp, ivtmp.64
	vbroadcastsd	(%rdx,%rbx), %ymm12	# MEM[(const struct complex &)_95], _153
	vbroadcastsd	(%rdx,%rcx), %ymm14	# MEM[(const struct complex &)_237], _151
	vbroadcastsd	8(%rdx,%rcx), %ymm11	# MEM[(const struct complex &)_237 + 8], _154
	vbroadcastsd	8(%rdx,%rbx), %ymm9	# MEM[(const struct complex &)_95 + 8], _156
	vbroadcastsd	(%rdx,%r10), %ymm13	# MEM[(const struct complex &)_82], _152
	vbroadcastsd	8(%rdx,%r10), %ymm10	# MEM[(const struct complex &)_82 + 8], _155
# fft_batch_real.cpp:144: 			for (int ihi = 0; ihi < NHI; ihi++) {
	xorl	%r10d, %r10d	# ihi
	.p2align 4
	.p2align 3
.L22:
	leaq	0(,%r9,8), %rcx	#, _641
	leaq	0(%r13,%r9), %r8	#, tmp411
	leaq	(%r15,%rcx), %rdx	#, ivtmp.39
	leaq	(%r14,%r8,8), %r8	#, _658
	addq	%rax, %rcx	# ivtmp.72, ivtmp.40
	.p2align 4
	.p2align 3
.L23:
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:314:   return (__m256d) ((__v4df)__A * (__v4df)__B);
	vmulpd	(%rcx,%rsi), %ymm10, %ymm1	# MEM[base: _296, index: _502, offset: 0B], _155, tmp422
	vmulpd	(%rcx,%rsi,2), %ymm11, %ymm4	# MEM[base: _296, index: _502, step: 2, offset: 0B], _154, tmp418
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:362:   return (__m256d) ((__v4df)__A - (__v4df)__B);
	vmovapd	(%rdx), %ymm5	# MEM[base: _286, offset: 0B], tmp667
# /usr/lib/gcc/x86_64-linux-gnu/9/include/fmaintrin.h:97:   return (__m256d)__builtin_ia32_vfmsubpd256 ((__v4df)__A, (__v4df)__B,
	vfmsub231pd	(%rdx,%rsi), %ymm13, %ymm1	# MEM[base: _286, index: _502, offset: 0B], _152, tmp421
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:314:   return (__m256d) ((__v4df)__A * (__v4df)__B);
	vmulpd	(%rcx,%rdi), %ymm9, %ymm0	# MEM[base: _296, index: _677, offset: 0B], _156, tmp426
# /usr/lib/gcc/x86_64-linux-gnu/9/include/fmaintrin.h:97:   return (__m256d)__builtin_ia32_vfmsubpd256 ((__v4df)__A, (__v4df)__B,
	vfmsub231pd	(%rdx,%rsi,2), %ymm14, %ymm4	# MEM[base: _286, index: _502, step: 2, offset: 0B], _151, tmp417
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:314:   return (__m256d) ((__v4df)__A * (__v4df)__B);
	vmulpd	(%rcx,%rsi,2), %ymm14, %ymm7	# MEM[base: _296, index: _502, step: 2, offset: 0B], _151, tmp416
# /usr/lib/gcc/x86_64-linux-gnu/9/include/fmaintrin.h:97:   return (__m256d)__builtin_ia32_vfmsubpd256 ((__v4df)__A, (__v4df)__B,
	vfmsub231pd	(%rdx,%rdi), %ymm12, %ymm0	# MEM[base: _286, index: _677, offset: 0B], _153, tmp425
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:314:   return (__m256d) ((__v4df)__A * (__v4df)__B);
	vmulpd	(%rcx,%rsi), %ymm13, %ymm3	# MEM[base: _296, index: _502, offset: 0B], _152, tmp420
	vmulpd	(%rcx,%rdi), %ymm12, %ymm2	# MEM[base: _296, index: _677, offset: 0B], _153, tmp424
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:143:   return (__m256d) ((__v4df)__A + (__v4df)__B);
	vaddpd	(%rdx), %ymm1, %ymm6	# MEM[base: _286, offset: 0B], tmp421, _97
# /usr/lib/gcc/x86_64-linux-gnu/9/include/fmaintrin.h:49:   return (__m256d)__builtin_ia32_vfmaddpd256 ((__v4df)__A, (__v4df)__B,
	vfmadd231pd	(%rdx,%rsi,2), %ymm11, %ymm7	# MEM[base: _286, index: _502, step: 2, offset: 0B], _154, tmp415
	vfmadd231pd	(%rdx,%rsi), %ymm10, %ymm3	# MEM[base: _286, index: _502, offset: 0B], _155, tmp419
	vfmadd231pd	(%rdx,%rdi), %ymm9, %ymm2	# MEM[base: _286, index: _677, offset: 0B], _156, tmp423
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:362:   return (__m256d) ((__v4df)__A - (__v4df)__B);
	vsubpd	%ymm1, %ymm5, %ymm1	# tmp421, tmp667, _99
	vmovapd	(%rcx), %ymm5	# MEM[base: _296, offset: 0B], tmp668
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:143:   return (__m256d) ((__v4df)__A + (__v4df)__B);
	vaddpd	(%rcx), %ymm3, %ymm8	# MEM[base: _296, offset: 0B], tmp419, _98
	vaddpd	%ymm4, %ymm0, %ymm15	# tmp417, tmp425, _101
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:362:   return (__m256d) ((__v4df)__A - (__v4df)__B);
	vsubpd	%ymm4, %ymm0, %ymm0	# tmp417, tmp425, _103
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:143:   return (__m256d) ((__v4df)__A + (__v4df)__B);
	vaddpd	%ymm15, %ymm6, %ymm4	# _101, _97, tmp427
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:362:   return (__m256d) ((__v4df)__A - (__v4df)__B);
	vsubpd	%ymm15, %ymm6, %ymm6	# _101, _97, tmp431
	vsubpd	%ymm3, %ymm5, %ymm3	# tmp419, tmp668, _100
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:143:   return (__m256d) ((__v4df)__A + (__v4df)__B);
	vaddpd	%ymm7, %ymm2, %ymm5	# tmp415, tmp423, _102
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:868:   *(__m256d *)__P = __A;
	vmovapd	%ymm4, (%rdx)	# tmp427, MEM[base: _286, offset: 0B]
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:362:   return (__m256d) ((__v4df)__A - (__v4df)__B);
	vsubpd	%ymm2, %ymm7, %ymm2	# tmp423, tmp415, _104
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:143:   return (__m256d) ((__v4df)__A + (__v4df)__B);
	vaddpd	%ymm5, %ymm8, %ymm4	# _102, _98, tmp428
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:362:   return (__m256d) ((__v4df)__A - (__v4df)__B);
	vsubpd	%ymm8, %ymm5, %ymm5	# _98, _102, tmp432
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:868:   *(__m256d *)__P = __A;
	vmovapd	%ymm4, (%rcx,%rdi)	# tmp428, MEM[base: _296, index: _677, offset: 0B]
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:143:   return (__m256d) ((__v4df)__A + (__v4df)__B);
	vaddpd	%ymm2, %ymm1, %ymm4	# _104, _99, tmp429
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:362:   return (__m256d) ((__v4df)__A - (__v4df)__B);
	vsubpd	%ymm2, %ymm1, %ymm1	# _104, _99, tmp433
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:868:   *(__m256d *)__P = __A;
	vmovapd	%ymm4, (%rdx,%rsi)	# tmp429, MEM[base: _286, index: _502, offset: 0B]
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:143:   return (__m256d) ((__v4df)__A + (__v4df)__B);
	vaddpd	%ymm0, %ymm3, %ymm4	# _103, _100, tmp430
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:362:   return (__m256d) ((__v4df)__A - (__v4df)__B);
	vsubpd	%ymm3, %ymm0, %ymm0	# _100, _103, tmp434
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:868:   *(__m256d *)__P = __A;
	vmovapd	%ymm4, (%rcx,%rsi,2)	# tmp430, MEM[base: _296, index: _502, step: 2, offset: 0B]
	vmovapd	%ymm6, (%rcx,%rsi)	# tmp431, MEM[base: _296, index: _502, offset: 0B]
	vmovapd	%ymm5, (%rdx,%rsi,2)	# tmp432, MEM[base: _286, index: _502, step: 2, offset: 0B]
	addq	$32, %rdx	#, ivtmp.39
	vmovapd	%ymm1, (%rcx)	# tmp433, MEM[base: _296, offset: 0B]
	vmovapd	%ymm0, -32(%rdx,%rdi)	# tmp434, MEM[base: _286, index: _677, offset: 0B]
	addq	$32, %rcx	#, ivtmp.40
# fft_batch_real.cpp:145: 				for (int l = 0; l < L; l += SIMD_SIZE) {
	cmpq	%rdx, %r8	# ivtmp.39, _658
	jne	.L23	#,
# fft_batch_real.cpp:144: 			for (int ihi = 0; ihi < NHI; ihi++) {
	incl	%r10d	# ihi
	addq	%r11, %r9	# _653, ivtmp.53
# fft_batch_real.cpp:144: 			for (int ihi = 0; ihi < NHI; ihi++) {
	cmpl	%r10d, %r12d	# ihi, NHI
	jne	.L22	#,
	decl	-56(%rbp)	# %sfp
	movq	-136(%rbp), %rdx	# %sfp, _613
	movl	-56(%rbp), %ecx	# %sfp, ivtmp.69
	addq	%rdx, -88(%rbp)	# _613, %sfp
	movq	-128(%rbp), %rdx	# %sfp, _606
	addq	%rdx, -80(%rbp)	# _606, %sfp
	movq	-104(%rbp), %rdx	# %sfp, _561
	addq	%rdx, -72(%rbp)	# _561, %sfp
# fft_batch_real.cpp:133: 		for (int k2 = 1; k2 < N2 / 2; k2++) {
	movl	-60(%rbp), %edx	# %sfp, k2
	addq	-120(%rbp), %rbx	# %sfp, ivtmp.66
	addq	-112(%rbp), %rax	# %sfp, ivtmp.72
	subl	%ecx, %edx	# ivtmp.69, k2
	cmpl	-64(%rbp), %edx	# %sfp, k2
	jl	.L25	#,
	jmp	.L10	#
.L11:
# fft_batch_real.cpp:105: 		if (N2 > 1) {
	cmpl	$1, -60(%rbp)	#, %sfp
	je	.L10	#,
	jmp	.L20	#
.L36:
	vzeroupper
# fft_batch_real.cpp:209: }
	subq	$-128, %rsp	#,
	popq	%rbx	#
	popq	%r10	#
	.cfi_remember_state
	.cfi_def_cfa 10, 0
	popq	%r12	#
	popq	%r13	#
	popq	%r14	#
	popq	%r15	#
	popq	%rbp	#
	leaq	-8(%r10), %rsp	#,
	.cfi_def_cfa 7, 8
	ret	
.L2:
	.cfi_restore_state
# fft_batch_real.cpp:7: 	int NHI = N / (N1 * N2);
	movl	-144(%rbp), %ebx	# %sfp, N
	testl	%ebx, %ebx	# N
	leal	7(%rbx), %eax	#, tmp323
# fft_batch_real.cpp:8: 	const auto& W = twiddles(N);
	movl	%ebx, %edi	# N,
# fft_batch_real.cpp:7: 	int NHI = N / (N1 * N2);
	cmovns	%ebx, %eax	# tmp323,, N, N
	sarl	$3, %eax	#, N
	movl	%eax, -60(%rbp)	# N, %sfp
# fft_batch_real.cpp:8: 	const auto& W = twiddles(N);
	call	_Z8twiddlesi@PLT	#
	movq	%rax, -96(%rbp)	# tmp454, %sfp
# fft_batch_real.cpp:20: 		for (int ihi = 0; ihi < NHI; ihi++) {
	cmpl	$7, %ebx	#, N
	jle	.L7	#,
# fft_batch_real.cpp:25: 				xi1 = xi0 + D;
	movslq	-140(%rbp), %rsi	# %sfp,
	movq	%rsi, %rax	#,
# fft_batch_real.cpp:25: 				xi1 = xi0 + D;
	leaq	0(,%rsi,8), %rdx	#, _9
	testl	%esi, %esi	# L
	jle	.L7	#,
	vmovapd	.LC0(%rip), %ymm7	#, tmp440
	vmovapd	.LC1(%rip), %ymm6	#, tmp441
	decl	%eax	# tmp330
	leal	0(,%rsi,8), %r12d	#, tmp328
	vmovapd	.LC2(%rip), %ymm8	#, tmp442
	movq	%rsi, %r9	# _8, tmp329
	leaq	(%rsi,%rsi,2), %r8	#, tmp335
	leaq	(%rsi,%rsi,4), %rdi	#, tmp339
	shrl	$2, %eax	#,
	imulq	$56, %rsi, %rsi	#, _8, _390
	movslq	%r12d, %r12	# tmp328, _373
	salq	$4, %r9	#, tmp329
	salq	$2, %rax	#, _382
	leaq	0(,%r12,8), %r14	#, _368
	addq	%r15, %r9	# x, ivtmp.131
	salq	$4, %r8	#, tmp336
	movq	%rax, -56(%rbp)	# _382, %sfp
	salq	$3, %rdi	#, tmp340
	xorl	%r10d, %r10d	# ivtmp.129
# fft_batch_real.cpp:20: 		for (int ihi = 0; ihi < NHI; ihi++) {
	xorl	%ebx, %ebx	# ihi
	leaq	32(%r15), %r13	#, tmp443
.L8:
	movq	-56(%rbp), %rcx	# %sfp, _382
	leaq	(%r15,%r10,8), %rax	#, ivtmp.116
	addq	%r10, %rcx	# ivtmp.129, tmp342
	leaq	0(%r13,%rcx,8), %r11	#, _377
	movq	%r9, %rcx	# ivtmp.131, ivtmp.118
.L6:
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:143:   return (__m256d) ((__v4df)__A + (__v4df)__B);
	vmovapd	(%rax), %ymm3	# MEM[base: _409, offset: 0B], tmp601
	addq	$32, %rcx	#, ivtmp.118
	vaddpd	(%rax,%rdx), %ymm3, %ymm10	# MEM[base: _409, index: _9, offset: 0B], tmp601, _202
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:362:   return (__m256d) ((__v4df)__A - (__v4df)__B);
	vsubpd	(%rax,%rdx), %ymm3, %ymm1	# MEM[base: _409, index: _9, offset: 0B], tmp602, _203
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:143:   return (__m256d) ((__v4df)__A + (__v4df)__B);
	vmovapd	-32(%rcx), %ymm3	# MEM[base: _405, offset: 0B], tmp603
	vaddpd	-32(%rcx,%rdx), %ymm3, %ymm13	# MEM[base: _405, index: _9, offset: 0B], tmp603, _204
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:362:   return (__m256d) ((__v4df)__A - (__v4df)__B);
	vsubpd	-32(%rcx,%rdx), %ymm3, %ymm5	# MEM[base: _405, index: _9, offset: 0B], tmp604, _205
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:143:   return (__m256d) ((__v4df)__A + (__v4df)__B);
	vmovapd	(%rax,%rdx,4), %ymm3	# MEM[base: _409, index: _9, step: 4, offset: 0B], tmp605
	vaddpd	(%rax,%rdi), %ymm3, %ymm4	# MEM[base: _409, index: _396, offset: 0B], tmp605, _208
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:362:   return (__m256d) ((__v4df)__A - (__v4df)__B);
	vsubpd	(%rax,%rdi), %ymm3, %ymm0	# MEM[base: _409, index: _396, offset: 0B], tmp606, _209
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:143:   return (__m256d) ((__v4df)__A + (__v4df)__B);
	vmovapd	(%rax,%r8), %ymm3	# MEM[base: _409, index: _402, offset: 0B], tmp607
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:362:   return (__m256d) ((__v4df)__A - (__v4df)__B);
	vsubpd	(%rax,%rsi), %ymm3, %ymm2	# MEM[base: _409, index: _390, offset: 0B], tmp608, _211
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:143:   return (__m256d) ((__v4df)__A + (__v4df)__B);
	vaddpd	(%rax,%rsi), %ymm3, %ymm11	# MEM[base: _409, index: _390, offset: 0B], tmp607, _210
	vaddpd	%ymm13, %ymm10, %ymm9	# _204, _202, _206
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:362:   return (__m256d) ((__v4df)__A - (__v4df)__B);
	vsubpd	%ymm13, %ymm10, %ymm10	# _204, _202, tmp355
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:314:   return (__m256d) ((__v4df)__A * (__v4df)__B);
	vmulpd	%ymm7, %ymm2, %ymm3	# tmp440, _211, tmp347
	vmulpd	%ymm6, %ymm2, %ymm2	# tmp441, _211, tmp351
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:143:   return (__m256d) ((__v4df)__A + (__v4df)__B);
	vaddpd	%ymm11, %ymm4, %ymm12	# _210, _208, _212
# /usr/lib/gcc/x86_64-linux-gnu/9/include/fmaintrin.h:49:   return (__m256d)__builtin_ia32_vfmaddpd256 ((__v4df)__A, (__v4df)__B,
	vfmadd231pd	%ymm0, %ymm7, %ymm3	# _209, tmp440, tmp345
# /usr/lib/gcc/x86_64-linux-gnu/9/include/fmaintrin.h:97:   return (__m256d)__builtin_ia32_vfmsubpd256 ((__v4df)__A, (__v4df)__B,
	vfmsub132pd	%ymm6, %ymm2, %ymm0	# tmp441, tmp351, tmp349
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:143:   return (__m256d) ((__v4df)__A + (__v4df)__B);
	vaddpd	%ymm12, %ymm9, %ymm2	# _212, _206, tmp353
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:868:   *(__m256d *)__P = __A;
	vmovapd	%ymm2, (%rax)	# tmp353, MEM[base: _409, offset: 0B]
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:143:   return (__m256d) ((__v4df)__A + (__v4df)__B);
	vaddpd	%ymm0, %ymm1, %ymm2	# tmp349, _203, tmp354
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:362:   return (__m256d) ((__v4df)__A - (__v4df)__B);
	vsubpd	%ymm0, %ymm1, %ymm1	# tmp349, _203, tmp356
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:143:   return (__m256d) ((__v4df)__A + (__v4df)__B);
	vaddpd	%ymm3, %ymm5, %ymm0	# tmp345, _205, tmp358
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:868:   *(__m256d *)__P = __A;
	vmovapd	%ymm2, (%rax,%rdx)	# tmp354, MEM[base: _409, index: _9, offset: 0B]
	vmovapd	%ymm10, -32(%rcx)	# tmp355, MEM[base: _405, offset: 0B]
	vmovapd	%ymm1, -32(%rcx,%rdx)	# tmp356, MEM[base: _405, index: _9, offset: 0B]
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:362:   return (__m256d) ((__v4df)__A - (__v4df)__B);
	vsubpd	%ymm12, %ymm9, %ymm1	# _212, _206, tmp357
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:868:   *(__m256d *)__P = __A;
	vmovapd	%ymm1, (%rax,%rdx,4)	# tmp357, MEM[base: _409, index: _9, step: 4, offset: 0B]
	vmovapd	%ymm0, (%rax,%rdi)	# tmp358, MEM[base: _409, index: _396, offset: 0B]
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:362:   return (__m256d) ((__v4df)__A - (__v4df)__B);
	vsubpd	%ymm11, %ymm4, %ymm0	# _210, _208, tmp359
	addq	$32, %rax	#, ivtmp.116
# fft_batch_real.cpp:70: 				_mm256_store_pd(xi6, -R5);
	vxorpd	%ymm8, %ymm0, %ymm0	# tmp442, tmp359, tmp360
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:868:   *(__m256d *)__P = __A;
	vmovapd	%ymm0, -32(%rax,%r8)	# tmp360, MEM[base: _409, index: _402, offset: 0B]
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:362:   return (__m256d) ((__v4df)__A - (__v4df)__B);
	vsubpd	%ymm5, %ymm3, %ymm0	# _205, tmp345, tmp362
# /usr/lib/gcc/x86_64-linux-gnu/9/include/avxintrin.h:868:   *(__m256d *)__P = __A;
	vmovapd	%ymm0, -32(%rax,%rsi)	# tmp362, MEM[base: _409, index: _390, offset: 0B]
# fft_batch_real.cpp:21: 			for (int l = 0; l < L; l += SIMD_SIZE) {
	cmpq	%rax, %r11	# ivtmp.116, _377
	jne	.L6	#,
# fft_batch_real.cpp:20: 		for (int ihi = 0; ihi < NHI; ihi++) {
	incl	%ebx	# ihi
	addq	%r12, %r10	# _373, ivtmp.129
	addq	%r14, %r9	# _368, ivtmp.131
# fft_batch_real.cpp:20: 		for (int ihi = 0; ihi < NHI; ihi++) {
	cmpl	-60(%rbp), %ebx	# %sfp, ihi
	jl	.L8	#,
.L7:
# fft_batch_real.cpp:77: 		NHI = N / (N1 * N2);
	movl	-144(%rbp), %ebx	# %sfp, N
# fft_batch_real.cpp:75: 		N2 *= N1;
	movl	$8, -60(%rbp)	#, %sfp
# fft_batch_real.cpp:77: 		NHI = N / (N1 * N2);
	testl	%ebx, %ebx	# N
	leal	31(%rbx), %eax	#, tmp326
	cmovns	%ebx, %eax	# tmp326,, N, N
	sarl	$5, %eax	#, N
	movl	%eax, %r12d	# N, NHI
	jmp	.L3	#
	.cfi_endproc
.LFE9575:
	.size	_Z14fft_batch_realPdii, .-_Z14fft_batch_realPdii
	.section	.rodata.cst32,"aM",@progbits,32
	.align 32
.LC0:
	.long	1719614413
	.long	-1075404642
	.long	1719614413
	.long	-1075404642
	.long	1719614413
	.long	-1075404642
	.long	1719614413
	.long	-1075404642
	.align 32
.LC1:
	.long	1719614413
	.long	1072079006
	.long	1719614413
	.long	1072079006
	.long	1719614413
	.long	1072079006
	.long	1719614413
	.long	1072079006
	.align 32
.LC2:
	.long	0
	.long	-2147483648
	.long	0
	.long	-2147483648
	.long	0
	.long	-2147483648
	.long	0
	.long	-2147483648
	.ident	"GCC: (Ubuntu 9.4.0-1ubuntu1~20.04.1) 9.4.0"
	.section	.note.GNU-stack,"",@progbits
	.section	.note.gnu.property,"a"
	.align 8
	.long	 1f - 0f
	.long	 4f - 1f
	.long	 5
0:
	.string	 "GNU"
1:
	.align 8
	.long	 0xc0000002
	.long	 3f - 2f
2:
	.long	 0x3
3:
	.align 8
4:
