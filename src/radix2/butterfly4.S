#include      "ymm.h"

              .global        butterfly4


              .text

butterfly4:   vaddpd         er2, er0, tr0
              vaddpd         ei2, ei0, ti0
              vsubpd         er2, er0, tr2
              vsubpd         ei2, ei0, ti2
              vaddpd         er3, er1, tr1
              vaddpd         ei3, ei1, ti1
              vsubpd         er3, er1, tr3
              vsubpd         ei3, ei1, ti3
              vaddpd         tr1, tr0, er0
              vaddpd         ti1, ti0, ei0
              vaddpd         ti3, tr2, er1
              vsubpd         tr3, ti2, ei1
              vsubpd         tr1, tr0, er2
              vsubpd         ti1, ti0, ei2
              vsubpd         ti3, tr2, er3
              vaddpd         tr3, ti2, ei3
              ret
