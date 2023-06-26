

              .global        reverse_neg


              .text


reverse_neg:  not            %rdi
              mov            %rsi, %rdx
              dec            %rdx
              and            %rdx, %rdi
              bsf            %rsi, %rcx
              dec            %rcx
              mov            $1, %rdx
              shl            %rcx, %rdx
L501:         test           %rdx, %rdi
              jz             L10001
              xor            %rdx, %rdi
              shr            %rdx
              jmp            L501
L10001:       or             %rdx, %rdi
              mov            %rdi, %rax
              ret

