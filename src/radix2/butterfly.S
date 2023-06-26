#define       X              %r8
#define       Y              %r9
#define       W              %r10
#define       N              %r11
#define       STRIDE         %r12
#define       k2             %r13
#define       k2rev          %r14
#define       ilo            %r15
#define       UR0            %ymm0
#define       UR1            %ymm1
#define       UI0            %ymm2
#define       UI1            %ymm3
#define       TR0            %ymm4
#define       TR1            %ymm5
#define       TI0            %ymm6
#define       TI1            %ymm7
#define       COS            %ymm8
#define       SIN            %ymm9
#define       SGN            %ymm10
#define       N2             -8(%rbp)
#define       NLO            -16(%rbp)
#define       STACK_SIZE     $16


              .global        butterfly


              .text
butterfly:    push           %rbp
              mov            %rsp, %rbp
              sub            STACK_SIZE, %rsp
              push           %r12
              push           %r13
              push           %r14
              push           %r15
              push           %rbx
              mov            %r8, STRIDE
              mov            %rcx, N
              mov            %rdx, W
              mov            %rsi, Y
              mov            %rdi, X
              mov            $1, %rax
              mov            %rax, N2
              mov            N, %rax
              shr            %rax
              mov            %rax, NLO
// ilo + NLO * (n1 + N1 * k2)
L5000:        xor            k2, k2
              xor            k2rev, k2rev
L10000:       xor            ilo, ilo
L12500:       mov            NLO, %rax

              inc            ilo
              cmp            ilo, NLO
              jne            L12500
              inc            k2
              cmp            k2, N2
              je             L14500
              dec            k2
              mov            N2, %rax
              shr            %rax
L13000:       cmp            %rax, k2rev
              jl             L14000
              sub            %rax, k2rev
              shr            %rax
              jmp            L13000
L14000:       add            %rax, k2rev
              inc            k2
              jmp            L10000
L14500:       mov            N2, %rax
              shl            %rax
              cmp            %rax, N
              je             L15000
              mov            %rax, N2
              mov            NLO, %rax
              shr            %rax
              mov            %rax, NLO
              jmp            L5000
L15000:       pop            %rbp
              pop            %rbx
              pop            %r15
              pop            %r14
              pop            %r13
              pop            %r12
              mov            %rbp, %rsp
              pop            %rbp
              ret
/*

butterfly:    push           %r12
              push           %r13
              push           %r14
              push           %r15
              push           %rbx
              mov            %r8, STRIDE
              mov            %rcx, N
              mov            %rdx, W
              mov            %rsi, Y
              mov            %rdi, X
              mov            N, %rax
              shl            $6, %rax
              add            STACK_SIZE, %rax
              push           %rbp
              mov            %rsp, %rbp
              sub            %rax, %rsp
              mov            %rsp, %rsi
              xor            %rax, %rax
              xor            %rbx, %rbx
L10000:       mov            %rax, %rcx
              mov            %rbx, %rdx
              shl            $5, %rcx
              shl            $6, %rdx
              vmovapd        (X, %rcx), UR0
              vmovapd        (Y, %rcx), UI0
              vmovapd        UR0, (%rsi, %rdx)
              vmovapd        UI0, 32(%rsi, %rdx)
              inc            %rax
L14000:       cmp            %rax, N
              je             L15000
              mov            N, %rcx
              shr            %rcx
              cmp            %rcx, %rbx
              jl             L14500
              sub            %rcx, %rbx
              shr            %rcx
              jmp            L14000
L14500:       add            %rcx, %rbx
              jmp            L10000
L15000:       mov            N, %rax
              shr            %rax
              mov            %rax, NHI
              mov            $1, N2
L16000:       cmp            N2, N
              je             L20000
              xor            k2, k2
L17000:       mov            k2, %rax
              or             N2, %rax
              shl            $4, %rax
              vbroadcastsd   0(W, %rax), COS
              vbroadcastsd   8(W, %rax), SIN
              xor            ihi, ihi
L18000:       mov            ihi, %rax
              imul           N2, %rax
              shl            %rax
              add            k2, %rax
              mov            N2, %rbx
              shl            $6, %rax
              shl            $6, %rbx
              add            %rax, %rbx
              vmovapd        (%rsi, %rax), UR0
              vmovapd        (%rsi, %rax), UR1
              vmovapd        32(%rsi, %rbx), UI0
              vmovapd        32(%rsi, %rbx), UI1
              vmovapd        UR0, TR0
              vmovapd        UI0, TI0
              vfmadd231pd    SIN, UI1, TR0
              vfnmadd231pd   SIN, UR1, TI0
              vfnmadd132pd   COS, TR0, UR1
              vfnmadd132pd   COS, TI0, UI1
              vfmsub132pd    two, UR1, UR0
              vfmsub132pd    two, UI1, UI0
              vmovapd        UR0, (%rsi, %rax)
              vmovapd        UR1, (%rsi, %rax)
              vmovapd        UI0, 32(%rsi, %rbx)
              vmovapd        UI1, 32(%rsi, %rbx)
              inc            ihi
              cmp            ihi, NHI
              jne            L18000
              inc            k2
              cmp            k2, N2
              jne            L17000
              shl            N2
              mov            NHI, %rax
              shr            %rax
              mov            %rax, NHI
              jmp            L16000
              xor            %rax, %rax
              xor            %rbx, %rbx
              shr            N2
L18100:       mov            %rbx, %rdx
              cmp            %rbx, N2
              jle            L18125
              mov            X, %rdi
              vmovapd        none, SGN
              jmp            L18150
L18125:       mov            Y, %rdi
              sub            N2, %rdx
              vmovapd        one, SGN
L18150:       mov            %rax, %rcx
              shl            $5, %rcx
              shl            $5, %rdx
              vmulpd         SGN, UI0, UI0
              vmovapd        (%rsi, %rcx), UR0
              sub            N, %rcx
              inc            %rcx
              neg            %rcx
              vmovapd        32(%rsi, %rcx), UI0
              vmovapd        UR0, (%rdi, %rdx)
              vmovapd        UI0, 32(%rdi, %rdx)
              inc            %rax
L18200:       cmp            %rax, N
              je             L20000
              mov            N, %rcx
              shr            %rcx
              cmp            %rcx, %rbx
              jl             L18300
              sub            %rcx, %rbx
              shr            %rcx
              jmp            L18200
L18300:       add            %rcx, %rbx
              jmp            L18100
L20000:       mov            %rbp, %rsp
              pop            %rbp
              pop            %rbx
              pop            %r15
              pop            %r14
              pop            %r13
              pop            %r12
              ret
*/
