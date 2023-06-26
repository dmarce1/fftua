#define       N1             %r12
#define       W              %r13
#define       COS0           %xmm0
#define       SIN0           %xmm1
#define       dCOS           %xmm2
#define       dSIN           %xmm3
#define       COS1           %xmm4
#define       SIN1           %xmm5
#define       COS2           %xmm6
#define       SIN2           %xmm7

              .global        twiddles


              .text


twiddles:     push           %r12
              push           %r13
              mov            %rdi, W
              vmovq          %rsi, COS0
              vmovq          %rdx, SIN0
              mov            %rcx, N1
              vmovq          none, dCOS
              vmovq          zero, dSIN
              vmovq          COS0, COS1
              vmovq          SIN0, SIN1
              shr            %rdi
L10000:       mov            %rdi, %rax
              shl            $4, %rax
              vmovq          COS1, (W, %rax)
              vmovq          SIN1, 8(W, %rax)
              shr            %rdi
              cmp            $0, %rdi
              je             L20000
              vmulsd         SIN1, SIN1, COS2
              vmulsd         COS1, SIN1, SIN2
              vfmsub231sd    COS1, COS1, COS2
              vfmadd231sd    SIN1, COS1, SIN2
              vmovq          COS2, COS1
              vmovq          SIN2, SIN1
              jmp            L10000
L20000:       cmp            $2, N1
              je             L40000
              mov            $2, %rdi
L30000:       mov            $1, %rdx
L35000:       mov            %rdi, %rax
              or             %rdx, %rax
              shl            $4, %rax
              vmovq          -16(W, %rax), COS1
              vmovq          -8(W, %rax), SIN1
              vmulsd         dSIN, SIN1, COS2
              vmulsd         dCOS, SIN1, SIN2
              vfmsub231sd    dCOS, COS1, COS2
              vfmadd231sd    dSIN, COS1, SIN2
              vmovq          COS2, 0(W, %rax)
              vmovq          SIN2, 8(W, %rax)
              inc            %rdx
              cmp            %rdx, %rdi
              jne            L35000
              shl            %rdi
              cmp            %rdi, N1
              je             L40000
              vaddsd         one, dCOS, dCOS
              vmulsd         half, dCOS, dCOS
              vsubsd         one, dCOS, dSIN
              vmulsd         none, dSIN, dSIN
              sqrtsd         dSIN, dSIN
              sqrtsd         dCOS, dCOS
              vmulsd         none, dSIN, dSIN
              jmp            L30000
L40000:       pop            %r13
              pop            %r12
              ret




