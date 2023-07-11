#define       logN           %r12
#define       LOCK           %r13
#define       W              %r14
#define       STACK_SIZE     $64
#define       N2m1           -8(%rbp)
#define       N              -16(%rbp)
#define       twopioN        -24(%rbp)
#define       cos            -32(%rbp)
#define       sin            -40(%rbp)
#define       N2             -48(%rbp)
#define       ptr            -56(%rbp)
#define       N1             $4
#define       logN1          $2

              .global        twiddles_fwd_complex2


              .text

twiddles_fwd_complex2:
              push           %rbp
              mov            %rsp, %rbp
              sub            STACK_SIZE, %rsp
              push           %rbx
              push           %r12
              push           %r13
              push           %r14
              push           %r15
              mov            %rdi, N
              bsr            N, logN
              mov            N, %rax
              shr            logN1, %rax
              mov            %rax, N2
              dec            %rax
              mov            %rax, N2m1
              lea            (locks), %rax
              lea            (%rax, logN, 8), LOCK
spin_loop:    xor            %rax, %rax
              xchg           (LOCK), %rax
              cmp            $1, %rax
              je             lock_acquired
              pause
              jmp            spin_loop
lock_acquired:lea            (Ws), %rax
              lea            (%rax, logN, 8), %rax
              mov            (%rax), %rax
              cmp            $0, %rax
              jne            load_from_mem
              lea            ptr, %rdi
              imul           $3, N, %rdx
              mov            $32, %rsi
              call           posix_memalign
              mov            ptr, W
              lea            ptr, %rdi
              emms
              fldpi
              fldpi
              faddp
              fildl          N
              fdivrp
              fchs
              fstpl          twopioN
              xor            %rcx, %rcx
main_loop:    mov            %rcx, %rdx
              shl            %rdx
              mov            %rcx, cos
              fildl          cos
              fmull          twopioN
              fcos
              fstpl          (W, %rdx, 8)
              mov            %rcx, sin
              fildl          sin
              fmull          twopioN
              fsin
              fstpl          8(W, %rdx, 8)
              cmp            %rcx, N2m1
              je             load_complete
              inc            %rcx
              jmp            main_loop
load_complete:lea            (Ws), %rax
              lea            (%rax, logN, 8), %rax
              mov            W, (%rax)
              jmp            done
load_from_mem:lea            (Ws), %rax
              lea            (%rax, logN, 8), %rax
              mov            (%rax), W
done:         mov            W, %rax
              mov            $1, %rcx
              xchg           (LOCK), %rcx
              pop            %r15
              pop            %r14
              pop            %r13
              pop            %r12
              pop            %rbx
              mov            %rbp, %rsp
              pop            %rbp
              ret



              .data

              .align         8
Ws:           .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
locks:        .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
