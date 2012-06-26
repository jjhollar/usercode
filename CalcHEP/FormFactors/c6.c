/*******************************
*    CalcHEP version 2.5.1   *
*******************************/
#include<math.h>
#define Real double
#include"num_out.h"
#include"num_in.h"

extern double DkappA,LambA; 
extern double a0w,aCw;

extern void C6_ext(Real * C);
void C6_ext(Real * C)
{
Real* V=va_ext;
Real S[4];                                                                  
     
S[0]=V[10]*V[10];
S[0]*=S[0];
C[10]=+1536*S[0];
S[1]=V[10]*V[10];

C[9]=+S[1]*(S[1]*(64*aCw+448*a0w)+640);
S[2]=aCw*aCw;

C[8]=+S[1]*(S[1]*(a0w*(24*aCw+48*a0w)+4*S[2])+64*aCw+288*a0w)+320;
C[7]=+S[1]*(a0w*(16*aCw+32*a0w)+2*S[2])+48*aCw+160*a0w;
C[6]=+a0w*(8*aCw+16*a0w)+S[2];
C[5]=+S[1]*(160*aCw+64*a0w);
C[4]=+aCw*(4*S[1]*aCw+64)+64*a0w;
C[3]=+8*S[2];
C[2]=+4*S[2];
C[1]=+128*S[0];
S[3]=V[8]*V[8];
S[3]*=S[3];
C[0]=+S[3];
}
