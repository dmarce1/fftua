/*
 * asm.hpp
 *
 *  Created on: Jun 9, 2023
 *      Author: dmarce1
 */

#ifndef ASM_HPP_
#define ASM_HPP_



extern "C" {
void fft_swap( double* X, int NHI, int NA, int NMID, int NB, int NLO);
}

#endif /* ASM_HPP_ */
