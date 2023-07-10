#include      "ymm.h"

              .global        butterfly4


              .text

tw_butterfly4:vmovapd        er0, tr
              vmovapd        ei0, ti
              vfmadd231pd    sin2, ei2, tr
              vfnmadd231pd   sin2, er2, ti
              vfnmadd132pd   cos2, tr, er2
              vfnmadd132pd   cos2, ti, ei2
              vfmsub132pd    two, er2, er0
              vfmsub132pd    two, ei2, ei0
              vmovapd        er1, tr
              vmovapd        ei1, ti
              vfmadd231pd    sin2, ei3, tr
              vfnmadd231pd   sin2, er3, ti
              vfnmadd132pd   cos2, tr, er3
              vfnmadd132pd   cos2, ti, ei3
              vfmsub132pd    two, er3, er1
              vfmsub132pd    two, ei3, ei1
              vmovapd        er0, tr
              vmovapd        ei0, ti
              vfmadd231pd    sin1, ei1, tr
              vfnmadd231pd   sin1, er1, ti
              vfnmadd132pd   cos1, tr, er1
              vfnmadd132pd   cos1, ti, ei1
              vfmsub132pd    two, er1, er0
              vfmsub132pd    two, ei1, ei0
              vmovapd        er2, tr
              vmovapd        ei2, ti
              vfmadd231pd    cos1, ei3, tr
              vfnmadd231pd   cos1, er3, ti
              vfmadd132pd    sin1, tr, er3
              vfmadd132pd    sin1, ti, ei3
              vfmsub132pd    two, er3, er2
              vfmsub132pd    two, ei3, ei2
              ret
