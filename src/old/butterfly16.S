#define       X              %rdi
#define       Y              %rsi
#define       C              %rdx
#define       S              %rcx
#define       er0            %ymm0
#define       er1            %ymm1
#define       er2            %ymm2
#define       er3            %ymm3
#define       ei0            %ymm4
#define       ei1            %ymm5
#define       ei2            %ymm6
#define       ei3            %ymm7
#define       tr             %ymm8
#define       ti             %ymm9
#define       cos1           %ymm10
#define       cos2           %ymm11
#define       sin1           %ymm12
#define       sin2           %ymm13

              .global        butterfly16


              .data

              .align         32
one:          .double        1.0;
zero:         .double        0.0;
two:          .double        2.0;
              .double        2.0;
              .double        2.0;
              .double        2.0;

              .text

butterfly16:  vmovq          (C), %xmm0
              vmovq          (S), %xmm1
              vmulsd         %xmm1, %xmm1, %xmm2
              vmulsd         %xmm0, %xmm1, %xmm3
              vfmsub231sd    %xmm0, %xmm0, %xmm2
              vfmadd231sd    %xmm1, %xmm0, %xmm3
              vmulsd         %xmm4, %xmm3, %xmm4
              vmulsd         %xmm3, %xmm3, %xmm5
              vfmsub231sd    %xmm3, %xmm2, %xmm4
              vfmadd231sd    %xmm4, %xmm2, %xmm5
              vmulsd         %xmm6, %xmm5, %xmm8
              vmulsd         %xmm5, %xmm5, %xmm9
              vfmsub231sd    %xmm5, %xmm4, %xmm8
              vfmadd231sd    %xmm6, %xmm4, %xmm9
              vbroadcastsd   %xmm8, %ymm8
              vbroadcastsd   %xmm9, %ymm9
              movddup        %xmm0, %xmm0
              movddup        %xmm1, %xmm1
              movddup        %xmm2, %xmm2
              movddup        %xmm3, %xmm3
              movddup        %xmm4, %xmm4
              movddup        %xmm5, %xmm5
              vshufpd        $0, one, %xmm0, %xmm0
              vshufpd        $0, zero, %xmm1, %xmm1
              vshufpd        $0, %xmm2, %xmm4, %xmm2
              vshufpd        $0, %xmm3, %xmm5, %xmm3
              vinsertf128    $0, %xmm2, %ymm0, %ymm0
              vinsertf128    $0, %xmm3, %ymm1, %ymm1
              vmulpd         %ymm1, %ymm9, %ymm10
              vmulpd         %ymm0, %ymm9, %ymm11
              vfmsub231pd    %ymm0, %ymm8, %ymm10
              vfmadd231pd    %ymm1, %ymm8, %ymm11
              vmulpd         %ymm3, %ymm9, %ymm12
              vmulpd         %ymm2, %ymm9, %ymm13
              vfmsub231pd    %ymm2, %ymm8, %ymm12
              vfmadd231pd    %ymm3, %ymm8, %ymm13
              vmulpd         %ymm5, %ymm9, %ymm14
              vmulpd         %ymm4, %ymm9, %ymm15
              vfmsub231pd    %ymm4, %ymm8, %ymm14
              vfmadd231pd    %ymm5, %ymm8, %ymm15
              vmovapd        %ymm0, %ymm9
              vmovapd        %ymm1, %ymm10
              vmovapd        32(X), er1
              vmovapd        64(X), er2
              vmovapd        96(X), er3
              vmovapd        32(Y), ei1
              vmovapd        64(Y), ei2
              vmovapd        96(Y), ei3
              vmulpd         %ymm9, er1, er0
              vmulpd         %ymm10, er1, ei0
              vfmsub231pd    %ymm9, ei1, er0
              vfmadd231pd    %ymm11, ei1, ei0
              vmovapd        er0, er1
              vmovapd        ei0, ei1
              vmulpd         %ymm10, er1, er0
              vmulpd         %ymm11, er1, ei0
              vfmsub231pd    %ymm10, ei1, er0
              vfmadd231pd    %ymm11, ei1, ei0
              vmovapd        er0, er1
              vmovapd        ei0, ei1
              vmulpd         %ymm12, er2, er0
              vmulpd         %ymm13, er2, ei0
              vfmsub231pd    %ymm12, ei2, er0
              vfmadd231pd    %ymm13, ei2, ei0
              vmovapd        er0, er2
              vmovapd        ei0, ei2
              vmulpd         %ymm14, er3, er0
              vmulpd         %ymm15, er3, ei0
              vfmsub231pd    %ymm14, ei3, er0
              vfmadd231pd    %ymm15, ei3, ei0
              vmovapd        er0, er3
              vmovapd        ei0, ei3
              vmovapd        (X), er0
              vmovapd        (Y), ei0
              ret
