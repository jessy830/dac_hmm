/*
* Listing 3 -- bitio.h
*
* This header file contains the function prototypes needed to use
* the bitstream i/o routines.
*
*/


#include <intrin.h>  //zy
#include <math.h>     //zy


#define BUF_SIZE	32U
#define putone(buf, i)	_bittestandset((long*)buf+(i>>5), (i&31))	
//#define getbit(buf, i)	_bittest(buf+(i>>5), (i&31))
#define getbit(buf, i)  _bittest((long*)(buf)+((i)>>5), ((i)&31))//zy
bool input_bit(unsigned&);
void output_bit(unsigned&, unsigned);
void clearbuf(void);
