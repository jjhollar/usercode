/*******************************
*    CalcHEP version 2.5.1   *
*******************************/
#include<math.h>
#define Real double
#include"num_out.h"
#include"num_in.h"

extern double DkappA,LambA; 
extern double a0w,aCw;

extern void C3_ext(Real * C);
void C3_ext(Real * C)
{
Real* V=va_ext;
Real S[6];                                                                  
     
S[0]=V[10]*V[10];
S[0]*=S[0];S[0]*=S[0];
C[16]=+192*S[0];
S[1]=LambA*LambA;

S[2]=V[10]*V[10];

S[3]=V[10]*V[10];
S[3]*=S[3];S[3]*=V[10]*V[10];
C[15]=+S[3]*(S[2]*(8*aCw+48*a0w)+DkappA*(32*LambA+16*DkappA)+16*S[1]+288);
S[4]=V[10]*V[10];
S[4]*=S[4];
C[14]=+S[4]*(S[2]*(DkappA*(a0w*(16*LambA+8*DkappA+32)+aCw*(6*LambA+3*DkappA+16))+
 LambA*(LambA*(3*aCw+8*a0w)+16*aCw+32*a0w)+24*aCw+80*a0w)+DkappA*(32*LambA+
 64*DkappA+16)+16-80*LambA);
C[13]=+S[2]*(S[2]*(DkappA*(a0w*(24*LambA+4*DkappA+32)+aCw*(6*LambA-3*DkappA+12))+
 LambA*(LambA*(11*aCw+28*a0w)+16*aCw+48*a0w)+26*aCw+48*a0w)+DkappA*(96*
 LambA+40*DkappA+128)+LambA*(32*LambA+160)+136);
C[12]=+S[2]*(DkappA*(a0w*(16*LambA+8*DkappA+32)+aCw*(4*LambA+2*DkappA+8))+LambA*(
 LambA*(14*aCw+24*a0w)+8*aCw+32*a0w)+12*aCw+40*a0w)+DkappA*(-8*DkappA-16)-
 8;
C[11]=+288*S[3];
C[10]=+S[4]*(S[2]*(12*aCw+16*a0w)+DkappA*(320-32*DkappA)+LambA*(64*LambA+320)+
 304);
C[9]=+S[2]*(S[2]*(LambA*(a0w*(24*LambA+32*DkappA+48)+aCw*(8*(DkappA+1)+10*
 LambA))+DkappA*(aCw*(44-4*DkappA)+64*a0w)+50*aCw+64*a0w)+DkappA*(192*LambA+128*
 DkappA+304)+LambA*(32*LambA+256)+240);
C[8]=+S[2]*(DkappA*(a0w*(48*LambA+8*DkappA+64)+aCw*(28*LambA+2*DkappA+24))+LambA*(
 LambA*(34*aCw+40*a0w)+56*aCw+96*a0w)+40*aCw+88*a0w);
C[7]=+S[4]*(DkappA*(32*LambA+16*DkappA+320)+S[2]*(8*aCw+48*a0w)+LambA*(16*LambA+
 192)+304);
C[6]=+S[2]*(S[2]*(aCw*(8*LambA-32*DkappA-12)+a0w*(-32*DkappA-8))+DkappA*(-96*
 LambA-160*DkappA-256)+LambA*(-64*LambA-128)-112);
C[5]=+S[2]*(LambA*(a0w*(16*LambA+32*DkappA+96)+aCw*(32*LambA+48*DkappA+112))+
 DkappA*(aCw*(4*DkappA+40)+32*a0w)+60*aCw+56*a0w)+DkappA*(40*DkappA+80)+40;
C[4]=+S[2]*(64*(DkappA*(-DkappA-1)-S[1])+S[2]*(12*aCw+8*a0w)-32*LambA);
C[3]=+S[2]*(aCw*(DkappA*(32*LambA+8*DkappA+40)+LambA*(16*LambA+88)+48)+a0w*(32*
 LambA+8))+DkappA*(48*DkappA+96)+48;
C[2]=+S[2]*(aCw*(DkappA*(8*LambA+4*DkappA+16)+LambA*(4*LambA+24)+16))+DkappA*(16*
 DkappA+32)+16;
C[1]=+16*S[3];
S[5]=V[8]*V[8];
S[5]*=S[5];
C[0]=+S[5];
}
