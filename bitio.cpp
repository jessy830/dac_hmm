#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include "coder.h"
#include "bitio.h"

unsigned buf[BUF_SIZE+4];/* This is the i/o buffer    */
void output_bit(unsigned &p, unsigned bit){
	if(bit) putone(buf, p);
	p++;
}
bool input_bit(unsigned &p){
	bool bit =  getbit(buf, p);
	p++;
	return bit;
}
void clearbuf(void){
	memset(buf, 0, sizeof(buf));
}