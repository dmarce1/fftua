#include      "ymm.h"

              .global        transpose_lo


              .text

transpose_lo: vunpcklpd      er1, er0, %ymm8
              vunpckhpd      er1, er0, %ymm9
              vunpcklpd      er3, er2, %ymm10
              vunpckhpd      er3, er2, %ymm11
              vperm2f128     $0x20, %ymm10, %ymm8, er0
              vperm2f128     $0x20, %ymm11, %ymm9, er2
              vperm2f128     $0x31, %ymm10, %ymm8, er1
              vperm2f128     $0x31, %ymm11, %ymm9, er3
              vunpcklpd      ei1, ei0, %ymm8
              vunpckhpd      ei1, ei0, %ymm9
              vunpcklpd      ei3, ei2, %ymm10
              vunpckhpd      ei3, ei2, %ymm11
              vperm2f128     $0x20, %ymm10, %ymm8, ei0
              vperm2f128     $0x20, %ymm11, %ymm9, ei2
              vperm2f128     $0x31, %ymm10, %ymm8, ei1
              vperm2f128     $0x31, %ymm11, %ymm9, ei3
              ret
