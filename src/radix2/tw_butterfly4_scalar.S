#define       er0            %xmm0
#define       er1            %xmm1
#define       er2            %xmm2
#define       er3            %xmm3
#define       ei0            %xmm4
#define       ei1            %xmm5
#define       ei2            %xmm6
#define       ei3            %xmm7
#define       cos1           %xmm8
#define       sin1           %xmm9
#define       cos2           %xmm10
#define       sin2           %xmm11
#define       tr             %xmm12
#define       ti             %xmm13
#define       two            %xmm14

              .global        tw_butterfly4_scalar


              .text

tw_butterfly4_scalar:
              vmovq          er0, tr
              vmovq          ei0, ti
              vfmadd231sd    sin2, ei2, tr
              vfnmadd231sd   sin2, er2, ti
              vfnmadd132sd   cos2, tr, er2
              vfnmadd132sd   cos2, ti, ei2
              vfmsub132sd    two, er2, er0
              vfmsub132sd    two, ei2, ei0
              vmovq          er1, tr
              vmovq          ei1, ti
              vfmadd231sd    sin2, ei3, tr
              vfnmadd231sd   sin2, er3, ti
              vfnmadd132sd   cos2, tr, er3
              vfnmadd132sd   cos2, ti, ei3
              vfmsub132sd    two, er3, er1
              vfmsub132sd    two, ei3, ei1
              vmovq          er0, tr
              vmovq          ei0, ti
              vfmadd231sd    sin1, ei1, tr
              vfnmadd231sd   sin1, er1, ti
              vfnmadd132sd   cos1, tr, er1
              vfnmadd132sd   cos1, ti, ei1
              vfmsub132sd    two, er1, er0
              vfmsub132sd    two, ei1, ei0
              vmovq          er2, tr
              vmovq          ei2, ti
              vfmadd231sd    cos1, ei3, tr
              vfnmadd231sd   cos1, er3, ti
              vfmadd132sd    sin1, tr, er3
              vfmadd132sd    sin1, ti, ei3
              vfmsub132sd    two, er3, er2
              vfmsub132sd    two, ei3, ei2
              ret
