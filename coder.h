/*
* Listing 1 -- coder.h
*
* This header file contains the constants, declarations, and
* prototypes needed to use the arithmetic coding routines.  These
* declarations are for routines that need to interface with the
* arithmetic coding stuff in coder.c
*
*/

/*
* Function prototypes.
*/

typedef unsigned char byte;
//#define entropy(p)	(p)*(1-(p))==0.0 ? 0 : (-(p)*log2f(p) - (1-(p))*log2f(1-(p)))
#define entropy(p)	(p)*(1.0-(p))==0.0? 0.0 : (-(p)*log(p)/log(2.0) - (1.0-(p))*log(1.0-(p))/log(2.0))  //zy
void initLUT(float p);
void initLUTdac(float p, float k);

void init_AE(void);
void encode_symbol(unsigned&, byte, bool);
void flush_AE(unsigned&);

void init_AD(unsigned&, unsigned&, unsigned&, unsigned&);
byte get_symbol(unsigned, unsigned, unsigned, bool, double*,double*,float);
void remove_symbol(unsigned&, unsigned&, unsigned&, unsigned&, byte, bool);