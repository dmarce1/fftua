#define       NHIm1          %rbx
#define       I0             %rsi
#define       J0             %rdi
#define       I              %r8
#define       J              %r9
#define       L              %r11
#define       X              %r12
#define       NHI            %r13
#define       NLO            %r15

              .global        scramble_hi

              .text

scramble_hi:  push           %rbx
              push           %r12
              push           %r13
              push           %r14
              push           %r15
              mov            %rdi, X
              mov            %rsi, NHI
              mov            %rdx, NLO
              mov            NHI, NHIm1
              dec            NHIm1
              xor            I, I
              xor            J, J
L0:           cmp            I, J
              jg             L5
              mov            I, I0
              mov            J, J0
              imul           NLO, I0
              imul           NLO, J0
              mov            I0, %rax
              mov            J0, %rdx
              mov            I0, L
              add            NLO, L
L1:           vmovapd        (X, %rax, 8), %ymm0
              vmovapd        (X, %rdx, 8), %ymm1
              vmovapd        %ymm1, (X, %rax, 8)
              vmovapd        %ymm0, (X, %rdx, 8)
              add            $4, %rax
              add            $4, %rdx
              cmp            %rax, L
              jg             L1
L5:           mov            J, %r10
              not            %r10
              and            NHIm1, %r10
              bsr            %r10, %rcx
              mov            $1, %r10
              shl            %rcx, %r10
              dec            %r10
              mov            J, %rcx
              and            %r10, %rcx
              inc            %r10
              or             %r10, %rcx
              mov            %rcx, J
              inc            I
              cmp            I, NHIm1
              jg             L0
              pop            %r15
              pop            %r14
              pop            %r13
              pop            %r12
              pop            %rbx
              ret

