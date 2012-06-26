/*******************************
*    CalcHEP version 2.5.1   *
*******************************/
#include<math.h>
#define Real double
#include"num_out.h"
#include"num_in.h"

extern double DkappA,LambA; 
extern double a0w,aCw;

extern void C4_ext(Real * C);
void C4_ext(Real * C)
{
Real* V=va_ext;
Real S[10];                                                                 
     
S[0]=V[10]*V[10];
S[0]*=S[0];S[0]*=S[0];S[0]*=V[10]*V[10]*V[10]*V[10];
C[17]=+48*S[0];
S[1]=LambA*LambA;

S[2]=V[10]*V[10];
S[2]*=S[2];S[2]*=S[2];
C[16]=+S[2]*(DkappA*(DkappA*(DkappA*(4*LambA+DkappA+16)+LambA*(6*LambA+32)+80)+LambA*(
 LambA*(4*LambA+16)+32)+128)+S[1]*(S[1]+16)+64);
S[3]=V[10]*V[10];
S[3]*=S[3];S[3]*=S[3];S[3]*=V[10]*V[10];
C[15]=+96*S[3];
C[14]=+S[2]*(DkappA*(64*(LambA+DkappA)+160)+LambA*(16*LambA+80)+100);
S[4]=V[10]*V[10];
S[4]*=S[4];S[4]*=V[10]*V[10];
C[13]=+S[4]*(DkappA*(DkappA*(DkappA*(-4*DkappA-40)+LambA*(24*LambA+48)-104)+LambA*(
 LambA*(32*LambA+88)+48)-128)+S[1]*(12*S[1]+88)-64);
C[12]=+S[2]*(DkappA*(64*LambA+32*DkappA+160)+LambA*(32*LambA+192)+76);
C[11]=+S[4]*(DkappA*(DkappA*(DkappA*(8*LambA+2*DkappA-32)+LambA*(12*LambA-32)-80)+LambA*
 (8*S[1]-32)-72)+LambA*(LambA*(2*S[1]+64)+16)-32);
S[5]=V[10]*V[10];
S[5]*=S[5];
C[10]=+S[5]*(DkappA*(DkappA*(DkappA*(16*(1-LambA)+4*DkappA)+LambA*(8*LambA-32)+32)+
 LambA*(LambA*(80*(LambA+1))-32)+32)+S[1]*(52*S[1]+80)+16);
C[9]=+S[4]*(DkappA*(-64*DkappA-16)+LambA*(48*LambA+112)+44);
S[6]=DkappA*DkappA;

C[8]=+S[5]*(DkappA*(LambA*(DkappA*(16*(DkappA-1)+40*LambA)+LambA*(48*(LambA+1))-16)+
 DkappA*(4*(S[6]+1))+24)+LambA*(LambA*(20*S[1]+112)+32)+20);
S[7]=V[10]*V[10];

C[7]=+S[1]*(S[7]*(DkappA*(64*(LambA-1)-32*DkappA)+96*S[1]-64));
C[6]=+S[5]*(LambA*(DkappA*(DkappA*(36*LambA+24*DkappA+96)+LambA*(24*LambA+160)+224)+
 LambA*(LambA*(6*LambA+64)+224)+176)+DkappA*(DkappA*(6*S[6]+32)+72)+40);
C[5]=+S[7]*(DkappA*(LambA*(LambA*(64*(LambA-1))+32*(-DkappA-1))+8*DkappA+16)+S[1]*(
 64*(S[1]-1))+8);
C[4]=+S[7]*(LambA*(LambA*(16*(DkappA*(DkappA+1)+S[1]+1)+LambA*(32*DkappA+64))+DkappA*(
 16*(-DkappA-1)))+DkappA*(4*DkappA+8)+4);
S[8]=LambA*LambA;
S[8]*=S[8];
C[3]=+64*S[8];
C[2]=+16*S[8];
C[1]=+4*S[2];
S[9]=V[8]*V[8];
S[9]*=S[9];
C[0]=+S[9];
}
