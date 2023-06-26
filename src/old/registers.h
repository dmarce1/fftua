/*
 * registers.h
 *
 *  Created on: Jun 24, 2023
 *      Author: dmarce1
 */
#ifndef REGISTERS_H_
#define REGISTERS_H_


#define  ilo          %r8
#define  k2           %r9
#define  NLO          %r10
#define  N2           %r11
#define  N            %r12
#define  X            %r13
#define  C            %r14
#define  X0           %r15
#define  er0          %ymm0
#define  er1          %ymm1
#define  er2          %ymm2
#define  er3          %ymm3
#define  ei0          %ymm4
#define  ei1          %ymm5
#define  ei2          %ymm6
#define  ei3          %ymm7
#define  tr0          %ymm8
#define  tr1          %ymm9
#define  tr2          %ymm10
#define  tr3          %ymm11
#define  ti0          %ymm12
#define  ti1          %ymm13
#define  ti2          %ymm14
#define  ti3          %ymm15
#define  ur0          %xmm0
#define  ur1          %xmm1
#define  ur2          %xmm2
#define  ur3          %xmm3
#define  ui0          %xmm4
#define  ui1          %xmm5
#define  ui2          %xmm6
#define  ui3          %xmm7
#define  sr0          %xmm8
#define  sr1          %xmm9
#define  sr2          %xmm10
#define  sr3          %xmm11
#define  si0          %xmm12
#define  si1          %xmm13
#define  si2          %xmm14
#define  si3          %xmm15
#define  SIMD_SIZE    $4
#define  STACK_SIZE   $72
#define  N0           -8(%rbp)
#define  cos3         -16(%rbp)
#define  cos2         -24(%rbp)
#define  cos1         -32(%rbp)
#define  cos0         -40(%rbp)
#define  sin3         -48(%rbp)
#define  sin2         -56(%rbp)
#define  sin1         -64(%rbp)
#define  sin0         -72(%rbp)

#endif /* REGISTERS_H_ */
