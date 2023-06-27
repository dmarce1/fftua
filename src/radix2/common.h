/*
 * common.h
 *
 *  Created on: Jun 25, 2023
 *      Author: dmarce1
 */

#ifndef COMMON_H_
#define COMMON_H_


#define       N1             $4
#define       ilo            %r8
#define       k2             %r9
#define       NLO            %r10
#define       N2             %r11
#define       N0             %r12
#define       X              %r13
#define       C              %r14
#define       S              %r15
#define       er0            %ymm0
#define       er1            %ymm1
#define       er2            %ymm2
#define       er3            %ymm3
#define       ei0            %ymm4
#define       ei1            %ymm5
#define       ei2            %ymm6
#define       ei3            %ymm7
#define       tr0            %ymm8
#define       ti0            %ymm9
#define       tr1            %ymm10
#define       ti1            %ymm11
#define       tr2            %ymm12
#define       ti2            %ymm13
#define       tr3            %ymm14
#define       ti3            %ymm15
#define       ur0            %xmm0
#define       ui0            %xmm1
#define       ur1            %xmm2
#define       ui1            %xmm3
#define       ur2            %xmm4
#define       ui2            %xmm5
#define       ur3            %xmm6
#define       ui3            %xmm7
#define       sr0            %xmm8
#define       si0            %xmm9
#define       sr1            %xmm10
#define       sr2            %xmm11
#define       si1            %xmm12
#define       si2            %xmm13
#define       sr3            %xmm14
#define       si3            %xmm15
#define       prmt_cntrl     $27
#define       STACK_SIZE     $16
#define       X0             -8(%rbp)
#define       N              -16(%rbp)



#endif /* COMMON_H_ */
