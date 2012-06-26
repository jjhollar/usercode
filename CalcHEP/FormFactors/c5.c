/*******************************
*    CalcHEP version 2.5.1   *
*******************************/
#include<math.h>
#define Real double
#include"num_out.h"
#include"num_in.h"

extern double DkappA,LambA; 
extern double a0w,aCw;

extern void C5_ext(Real * C);
void C5_ext(Real * C)
{
Real* V=va_ext;
Real S[6];                                                                  
     
S[0]=V[10]*V[10];
S[0]*=S[0];S[0]*=S[0];
C[15]=+192*S[0];
S[1]=LambA*LambA;

S[2]=V[10]*V[10];

S[3]=V[10]*V[10];
S[3]*=S[3];S[3]*=V[10]*V[10];
C[14]=+S[3]*(S[2]*(8*aCw+48*a0w)+DkappA*(32*LambA+16*DkappA)+16*S[1]);
S[4]=V[10]*V[10];
S[4]*=S[4];
C[13]=+S[4]*(S[2]*(DkappA*(a0w*(16*LambA+8*DkappA+32)+aCw*(6*LambA+3*DkappA+16))+
 LambA*(LambA*(3*aCw+8*a0w)+16*aCw+32*a0w)+28*aCw+48*a0w)+16*(DkappA*(
 DkappA+1)+1)+LambA*(48*(LambA+1)));
C[12]=+S[2]*(S[2]*(DkappA*(DkappA*(aCw+4*a0w)+LambA*(-2*aCw-8*a0w))+S[1]*(
 aCw+4*a0w))+DkappA*(8*DkappA+16)+8);
C[11]=+288*S[3];
C[10]=+S[4]*(DkappA*(64*(LambA+DkappA)+320)+S[2]*(4*aCw+80*a0w)+LambA*(64-32*
 LambA)+304);
C[9]=+S[2]*(S[2]*(LambA*(a0w*(24*LambA+32*DkappA+48)+aCw*(10*LambA+8*DkappA+24))+
 aCw*(DkappA*(-4*DkappA-20)-10)+24*a0w)+LambA*(96*(LambA+1))+16*(1-DkappA));
C[8]=+S[2]*(DkappA*(DkappA*(2*aCw+8*a0w)+LambA*(-4*aCw-16*a0w))+S[1]*(2*aCw+
 8*a0w));
C[7]=+S[4]*(DkappA*(32*LambA+16*DkappA+320)+S[2]*(8*aCw+48*a0w)+LambA*(16*LambA+
 192)+304);
C[6]=+S[2]*(S[2]*(aCw*(8*LambA-32*DkappA-48)+a0w*(32*(-DkappA-1)))+DkappA*(32*
 DkappA-96*LambA-64)+LambA*(128*LambA-32)-112);
C[5]=+S[2]*(aCw*(DkappA*(4*DkappA+16)+LambA*(8*(LambA-1))+12)+a0w*(32*(DkappA*(
 LambA+1)+1)+16*S[1]))+DkappA*(-8*DkappA-16)-8;
C[4]=+S[2]*(64*(DkappA*(-DkappA-1)-S[1])+S[2]*(12*aCw+8*a0w)-32*LambA);
C[3]=+S[2]*(aCw*(DkappA*(8*DkappA+24)+8*LambA+16)+a0w*(-32*LambA-8))+DkappA*(16*
 DkappA+32)+16;
C[2]=+S[2]*(aCw*(DkappA*(8*LambA+4*DkappA+16)+LambA*(4*LambA+24)+16))+DkappA*(16*
 DkappA+32)+16;
C[1]=+16*S[3];
S[5]=V[8]*V[8];
S[5]*=S[5];
C[0]=+S[5];
}
