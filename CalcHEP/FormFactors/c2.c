/*******************************
*    CalcHEP version 2.5.1   *
*******************************/
#include<math.h>
#define Real double
#include"num_out.h"
#include"num_in.h"

extern double DkappA,LambA; 
extern double a0w,aCw;

extern void C2_ext(Real * C);
void C2_ext(Real * C)
{
Real* V=va_ext;
Real S[11];                                                                 
     
S[0]=V[10]*V[10];
S[0]*=S[0];S[0]*=S[0];S[0]*=V[10]*V[10]*V[10]*V[10];
C[19]=+48*S[0];
S[1]=V[10]*V[10];
S[1]*=S[1];S[1]*=S[1];S[1]*=V[10]*V[10];
C[18]=+48*S[1];
S[2]=V[10]*V[10];
S[2]*=S[2];S[2]*=S[2];
C[17]=+S[2]*(DkappA*(DkappA*(DkappA*(8*LambA+2*DkappA+24)+LambA*(12*LambA+56)+100)+
 LambA*(LambA*(8*LambA+40)+104)+128)+LambA*(LambA*(LambA*(2*LambA+8)+60)+40)+50);
S[3]=DkappA*DkappA;

S[4]=V[10]*V[10];
S[4]*=S[4];S[4]*=V[10]*V[10];
C[16]=+S[4]*(LambA*(LambA*(DkappA*(8*LambA+6*DkappA+24)+LambA*(3*LambA+8)+44)+DkappA*(
 16*DkappA+60)+52)+DkappA*(DkappA*(32-S[3])+68)+36);
S[5]=LambA*LambA;

S[6]=V[10]*V[10];
S[6]*=S[6];
C[15]=+S[6]*(S[5]*(4*(LambA*(LambA+1)+S[3])+12*DkappA+16)+DkappA*(2*DkappA+4)+2);
C[14]=+S[6]*(LambA*(DkappA*(DkappA*(-8*DkappA-40)+LambA*(16*(LambA-1))-48)+LambA*(
 LambA*(10*LambA-8)-8)-8)+DkappA*(DkappA*(2*S[3]-22)-32)-10);
S[7]=V[10]*V[10];

C[13]=+S[5]*(S[7]*(LambA*(16*LambA+8)-8*DkappA));
C[12]=+S[2]*(DkappA*(64*LambA+32*DkappA+160)+LambA*(32*LambA+192)+172);
C[11]=+S[4]*(LambA*(DkappA*(DkappA*(24*LambA+16*DkappA+48)+LambA*(16*LambA+64)+56)+
 LambA*(LambA*(4*LambA+16)+100)+56)+S[3]*(4*S[3]-20)+26);
C[10]=+S[6]*(LambA*(DkappA*(DkappA*(24*(-LambA-DkappA)-120)-144*LambA-240)+LambA*(
 LambA*(6*LambA-56)-184)-168)+DkappA*(DkappA*(-2*S[3]-22)-72)-50);
C[9]=+S[7]*(S[5]*(DkappA*(-8*LambA-4*DkappA-32)+LambA*(4*LambA+32)+16)+DkappA*(6*
 DkappA+12)+6);
S[8]=LambA*LambA;
S[8]*=S[8];
C[8]=+16*S[8];
S[9]=DkappA*DkappA;
S[9]*=DkappA;
C[7]=+S[6]*(LambA*(DkappA*(DkappA*(48*LambA+32*DkappA+160)+LambA*(32*LambA+256)+384)+
 LambA*(LambA*(8*LambA+96)+352)+320)+DkappA*(8*S[9]+80)+80);
C[6]=+S[7]*(S[5]*(DkappA*(16*LambA+8*DkappA+48)+LambA*(24*LambA-48)-32)+DkappA*(-12*
 DkappA-24)-12);
C[5]=+S[6]*(LambA*(DkappA*(DkappA*(24*LambA+16*DkappA+80)+LambA*(16*LambA+128)+192)+
 LambA*(LambA*(4*LambA+48)+176)+160)+DkappA*(4*S[9]+40)+40);
C[4]=+S[7]*(S[5]*(DkappA*(8*LambA+4*DkappA+24)+LambA*(12*LambA-24)-16)+DkappA*(-6*
 DkappA-12)-6);
C[3]=+24*S[8];
C[2]=+8*S[8];
C[1]=+2*S[2];
S[10]=V[8]*V[8];
S[10]*=S[10];
C[0]=+S[10];
}
