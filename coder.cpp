/*

*/

#include <stdio.h>
#include <float.h>
#include "coder.h"
#include "bitio.h"
#include <math.h>
#include <intrin.h>
#include <cstdlib>
#include <iostream>

using namespace std;

typedef unsigned char byte;

/*
* These four variables define the current state of the arithmetic coder/decoder.  
*/
#define ACBITS	16
#define RANGE	(1<<ACBITS) 
#define HRANGE	(RANGE>>1)  
#define QRANGE	(HRANGE>>1) 
unsigned low, high, code, underflow;  
unsigned thresh[HRANGE], thresh_dac[HRANGE][2];
extern double *pdf;
extern double p;

extern double ys[100001];


void initLUT(float p){
	for(unsigned i=0; i<HRANGE; i++)
		thresh[i] = (unsigned)((i+HRANGE+1) * (1-p));		
}
void initLUTdac(float p, float k){
	for(unsigned i=0; i<HRANGE; i++){	
		thresh_dac[i][0] = (unsigned)((i+HRANGE+1) * powf(1-p,k));	
		thresh_dac[i][1] = (unsigned)((i+HRANGE+1) * (1-powf(p,k)));
	}
}


/*
* This routine must be called to initialize the encoding process.
* The high register is initialized to all 1s, and it is assumed that
* it has an infinite string of 1s to be shifted into the lower bit
* positions when needed.
*/
void init_AE(void){
	low = underflow = 0;
	high = RANGE - 1;
	clearbuf();
}

/*
* This routine is called to encode a symbol.  The symbol is passed
* in the SYMBOL structure as a low count, a high count, and a range,
* instead of the more conventional probability ranges.  The encoding
* process takes two steps.  First, the values of high and low are
* updated to take into account the range restriction created by the
* new symbol.  Then, as many bits as possible are shifted out to
* the output stream.  Finally, high and low are stable again and
* the routine returns.
*/
void encode_symbol(unsigned &ptr, byte c, bool dac){
	bool overflow = _bittest((long*)&high, ACBITS);// test the ACBITS-th bit
	unsigned index = high - low - HRANGE;
	if(dac){       
		if(c){
			low += thresh_dac[index][1];
		}else high = low + thresh_dac[index][0] - 1;
	}else{//end of NTAILS(32)
		if(c){
			low += thresh[index];//high--;
		}else high = low + thresh[index] - 1;
	}

	if(_bittest((long*)&high, ACBITS) != _bittest((long*)&low, ACBITS)){// test the ACBITS-th bit, high: 1x...x, low: 0x...x 
		while((high-low)<HRANGE){// drop underflow bits, i.e., the (ACBITS+1)-th bit
			underflow++;
			low  -= HRANGE;	low  <<= 1;
			high -= HRANGE;	high <<= 1;	high |= 1;
		}
	}else{
		if(overflow){// switch to ACBITS
			output_bit(ptr, _bittest((long*)&high, ACBITS));// output the MSB, i.e., the (ACBITS+1)-th bit
			while(underflow){// output underflow bits
				output_bit(ptr, !_bittest((long*)&high, ACBITS));
				underflow--;
			}
		}

		while(_bittest((long*)&high, ACBITS-1) == _bittest((long*)&low, ACBITS-1)){/* If this test passes, it means that the MSDigits match, and can be sent to the output stream. */
			output_bit(ptr, _bittest((long*)&high, ACBITS-1));
			while(underflow){
				output_bit(ptr, !_bittest((long*)&high, ACBITS-1));
				underflow--;
			}
			low  <<= 1;
			high <<= 1;	high |= 1;
		}
		low  &= (RANGE-1);
		high &= (RANGE-1);

		while(!_bittest((long*)&high, ACBITS-2) && _bittest((long*)&low, ACBITS-2)){	/* If this test passes, the numbers are in danger of underflow, because the MSDigits don't match, and the 2nd digits are just one apart. */
			underflow++;
			low  -= QRANGE;	low  <<= 1;
			high -= QRANGE;	high <<= 1;	high |= 1;
		}

		if((high-low)<HRANGE){// switch to (ACBITS+1) 
			low  <<= 1;
			high <<= 1;	high |= 1;
		}
	}
}

/*
* At the end of the encoding process, there are still significant
* bits left in the high and low registers.  We output two bits,
* plus as many underflow bits as are necessary.
*/
void flush_AE(unsigned &ptr){
	if(_bittest((long*)&high, ACBITS)){/* test the ACBITS-th bit */
		output_bit(ptr, 1);// output the ACCBITS-th bit
		while(underflow){// output underflow bits
			output_bit(ptr, 0);
			underflow--;
		}
		for(unsigned i=0; i<ACBITS; i++){
			output_bit(ptr, _bittest((long*)&high, ACBITS-1));
			high <<= 1;
		}
	}else{/* ACBITS */
		output_bit(ptr, _bittest((long*)&high, ACBITS-1));	// output the (ACCBITS-1)-th bit
		while(underflow){// output underflow bits
			output_bit(ptr, !_bittest((long*)&high, ACBITS-1));
			underflow--;
		}
		for(unsigned i=1; i<ACBITS; i++){
			output_bit(ptr, _bittest((long*)&high, ACBITS-2));
			high <<=1;
		}
	}
}

/*
* This routine is called to initialize the state of the arithmetic
* decoder.  This involves initializing the high and low registers
* to their conventional starting values, plus reading the first
* 16 bits from the input stream into the code value.
*/
void init_AD(unsigned &low, unsigned &high, unsigned &code, unsigned &ptr){
	high = (RANGE-1);
	ptr = code = low = 0;
	for(unsigned i=0; i<ACBITS; i++){
		code <<= 1;
		code |= input_bit(ptr);   //bool bit =  getbit(buf, p)
	}
}

byte get_symbol(unsigned low, unsigned high, unsigned code, bool dac, double *weight0, double *weight1, float k) {
	if (dac){
		unsigned range = (high - low + 1) ;
		
		if(code < (low + thresh_dac[high-low-HRANGE][1]))		
			return 0;
		else if(code >= (low + thresh_dac[high-low-HRANGE][0]))	
			return 1;
		else {
			double x = (code - low)/double(range);
			double u0 = x/(powf((1.0 - p), k)), u1 = (x - (1.0 - (powf(p, k))))/(powf(p, k));
			*weight0 = ys[int(u0 * 100000)];
			*weight1 = ys[int(u1 * 100000)];
			*weight0 = *weight0 / (*weight0 + *weight1);
			*weight1 = 1.0 - *weight0;
			if (*weight0 < 0) { 
				*weight0 = 0.00000000001;
				*weight1 = 1 - *weight0;
			}
			if (*weight1 < 0) {
				*weight1 = 0.00000000001;
				*weight0 = 1 - *weight1;
			}
			
			return 'x';	// 'ambiguity'
		}
			
	}else{
		if(code < (low + thresh[high-low-HRANGE]))	return 0;
		else return 1;
	}
	
}

/*
* Just figuring out what the present symbol is doesn't remove
* it from the input bit stream.  After the character has been
* decoded, this routine has to be called to remove it from the
* input stream.
*/
void remove_symbol(unsigned &low, unsigned &high, unsigned &code, unsigned &ptr, byte c, bool dac){
	unsigned index = high - low - HRANGE;
	if(dac){
		if(c){// bit 1
			low += thresh_dac[index][1];//high--;
		}else high = low + thresh_dac[index][0] - 1;
	}else{
		if(c){// bit 1
			low += thresh[index];//high--;
		}else high = low + thresh[index] - 1;
	}

	if(_bittest((long*)&high, ACBITS) != _bittest((long*)&low, ACBITS)){
		while((high-low)<HRANGE){	// drop underflow bits
			low  -= HRANGE;	low  <<= 1;
			high -= HRANGE;	high <<= 1;	high |= 1;
			code -= HRANGE;	code <<= 1;	code |= input_bit(ptr);
		} 
	} else {
		while(_bittest((long*)&high, ACBITS-1) == _bittest((long*)&low, ACBITS-1)){/* If the MSDigits match, the bits will be shifted out. */
			low  <<= 1;
			high <<= 1;	high |= 1;
			code <<= 1;	code |= input_bit(ptr);
		}
		low  &= (RANGE-1);
		high &= (RANGE-1);
		code &= (RANGE-1);

		while(!_bittest((long*)&high, ACBITS-2) && _bittest((long*)&low, ACBITS-2)){/* If underflow is threatening, shift out the 2nd MSDigit. */
			low  -= QRANGE;	low  <<= 1;
			high -= QRANGE;	high <<= 1;	high |= 1;
			code -= QRANGE;	code <<= 1;	code |= input_bit(ptr);
		}

		if((high-low)<HRANGE){	/* switch to ACBITS+1 */
			low  <<= 1;
			high <<= 1;	high |= 1;
			code <<= 1;	code |= input_bit(ptr);
		}
	}
}