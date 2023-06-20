#define  X            %r8
#define  M            %r9
#define  Ix           %r10
#define  Iy           %r11
#define  N1N2         %r15
#define  M0           $4
#define  STACK_SIZE   $16
#define  N1           -8(%rbp)
#define  N2           -16(%rbp)

         .global      fft_transpose_recursive

         .data

         .text


fft_transpose_recursive:
         enter        STACK_SIZE, $0
         push         %r12
         push         %r13
         push         %r14
         push         %r15
         push         %rbx
         mov          %rdi, X
         mov          %rsi, N1
         mov          %rdx, N2
         imul         N1, %rdx
         mov          %rdx, N1N2
         xor          Ix, Ix
         xor          Iy, Iy
         mov          N1, M
         call         recurse
         pop          %rbx
         pop          %r15
         pop          %r14
         pop          %r13
         pop          %r12
         leave
         ret
recurse: cmp          M0, M
         jg           further
         xor          %r12, %r12
loop0:   mov          Ix, %rsi
         mov          Iy, %rdi
         imul         N2, %rsi
         imul         N2, %rdi
         add          %r12, %rsi
         add          %r12, %rdi
         imul         N1, %rsi
         imul         N1, %rdi
         add          Iy, %rsi
         add          Ix, %rdi
         xor          %r13, %r13
loop1:   xor          %r14, %r14
loop2:   mov          %r13, %rax
         mov          %r14, %rbx
         imul         N1N2, %rax
         imul         N1N2, %rbx
         add          %r14, %rax
         add          %r13, %rbx
         add          %rsi, %rax
         add          %rdi, %rbx
         cmp          %rax, %rbx
         jle          skip
         lea          (X, %rax, 8), %rax
         lea          (X, %rbx, 8), %rbx
         mov          (%rax), %rcx
         mov          (%rbx), %rdx
         mov          %rcx, (%rbx)
         mov          %rdx, (%rax)
         mov          %r13, %rax
         mov          %r14, %rbx
skip:    inc          %r14
         cmp          M0, %r14
         jl           loop2
endloop: inc          %r13
         cmp          M0, %r13
         jl           loop1
         inc          %r12
         cmp          %r12, N2
         jg           loop0
         ret
further: cmp          Ix, Iy
         je           xyeq
         shr          M
         call         recurse
         add          M, Iy
         call         recurse
         sub          M, Iy
         add          M, Ix
         call         recurse
         add          M, Iy
         call         recurse
         sub          M, Iy
         sub          M, Ix
         shl          M
         ret
xyeq:    shr          M
         call         recurse
         add          M, Iy
         call         recurse
         add          M, Ix
         call         recurse
         sub          M, Iy
         sub          M, Ix
         shl          M
         ret






