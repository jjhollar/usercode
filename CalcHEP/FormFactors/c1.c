/*******************************
*    CalcHEP version 2.5.1   *
*******************************/
#include<math.h>
#define Real double
#include"num_out.h"
#include"num_in.h"
//#include <stdio.h> // needed for printf tests

extern double DkappA,LambA; 
extern double a0w,aCw;

extern void C1_ext(Real * C);
void C1_ext(Real * C)
{
 //print test
 //printf("in c1.c: DkappA=%g, LambdA=%g, a0w=%g, aCw=%g\n",DkappA,LambA,a0w,aCw);
 //
Real* V=va_ext;
Real S[10];                                                                 
     
S[0]=V[10]*V[10];
S[0]*=S[0];S[0]*=S[0];S[0]*=V[10]*V[10]*V[10]*V[10];
C[24]=+48*S[0];
S[1]=LambA*LambA;

S[2]=V[10]*V[10];
S[2]*=S[2];S[2]*=S[2];
C[23]=+S[2]*(DkappA*(DkappA*(DkappA*(4*LambA+DkappA+16)+LambA*(6*LambA+32)-16)+LambA*(
 LambA*(4*LambA+16)-96)-192)+LambA*(LambA*(S[1]-32)-272)-112);
S[3]=V[10]*V[10];
S[3]*=S[3];S[3]*=V[10]*V[10];
C[22]=+S[3]*(DkappA*(DkappA*(DkappA*(8*(-LambA-1)-6*DkappA)+LambA*(12*LambA+80)+40)+
 LambA*(LambA*(24*LambA+88)+80)-40)+LambA*(LambA*(10*S[1]-24)-128)-76);
S[4]=V[10]*V[10];
S[4]*=S[4];
C[21]=+S[4]*(DkappA*(DkappA*(DkappA*(6*DkappA-8*LambA+16)+LambA*(4*LambA+80)+60)+LambA*(
 LambA*(56*LambA+192)+208)+80)+LambA*(LambA*(LambA*(38*LambA+64)+192)+144)+36);
S[5]=V[10]*V[10];

C[20]=+S[5]*(LambA*(LambA*(DkappA*(16*(1-DkappA)+32*LambA)+LambA*(48*LambA+64)+16)+
 DkappA*(16*(DkappA+1)))+DkappA*(-4*DkappA-8)-4);
S[6]=V[10]*V[10];
S[6]*=S[6];S[6]*=S[6];S[6]*=V[10]*V[10];
C[19]=+96*S[6];
C[18]=+S[2]*(DkappA*(192*LambA+128*DkappA+480)+LambA*(80*LambA+464)+252);
C[17]=+S[3]*(DkappA*(DkappA*(DkappA*(24-16*LambA-8*DkappA)+112*LambA+248)+LambA*(LambA*(
 16*LambA+88)+112)+64)+LambA*(LambA*(8*S[1]-184)-368)-132);
C[16]=+S[4]*(DkappA*(DkappA*(DkappA*(16*LambA+20*DkappA+32)+LambA*(40*LambA+368)+180)+
 LambA*(LambA*(112*LambA+656)+880)+280)+LambA*(LambA*(LambA*(68*LambA+256)+720)+
 608)+132);
C[15]=+S[5]*(LambA*(LambA*(DkappA*(96*LambA-16*DkappA+144)+LambA*(112*LambA+320)+
 144)+DkappA*(48*(DkappA+1)))+DkappA*(-12*DkappA-24)-12);
C[14]=+S[2]*(DkappA*(64*LambA+32*DkappA+160)+LambA*(32*LambA+192)+76);
C[13]=+S[3]*(DkappA*(DkappA*(DkappA*(8*LambA+2*DkappA-32)+LambA*(12*LambA-32)-272)+
 LambA*(8*S[1]-32)-120)+LambA*(LambA*(2*S[1]+208)+352)+100);
C[12]=+S[4]*(DkappA*(DkappA*(DkappA*(80*LambA+28*DkappA+16)+LambA*(104*LambA+592)+212)+
 LambA*(LambA*(80*LambA+896)+1360)+392)+LambA*(LambA*(LambA*(28*LambA+384)+1088)+
 960)+196);
C[11]=+S[5]*(LambA*(LambA*(DkappA*(128*LambA+64*DkappA+352)+LambA*(64*LambA+640)+
 352)+DkappA*(32*(DkappA+1)))+DkappA*(-8*DkappA-16)-8);
C[10]=+S[3]*(DkappA*(-64*DkappA-16)+LambA*(48*LambA+112)+44);
S[7]=DkappA*DkappA;

C[9]=+S[4]*(LambA*(DkappA*(DkappA*(104*LambA+80*DkappA+400)+LambA*(48*LambA+592)+912)+
 LambA*(LambA*(4*LambA+256)+784)+672)+DkappA*(DkappA*(20*S[7]+124)+264)+140);
C[8]=+S[5]*(DkappA*(LambA*(LambA*(128*(LambA+DkappA)+352)+32*(-DkappA-1))+8*DkappA+16)+
 S[1]*(640*LambA+352)+8);
S[8]=LambA*LambA;
S[8]*=S[8];
C[7]=+64*S[8];
C[6]=+S[4]*(LambA*(DkappA*(DkappA*(36*LambA+24*DkappA+96)+LambA*(24*LambA+160)+224)+
 LambA*(LambA*(6*LambA+64)+224)+176)+DkappA*(DkappA*(6*S[7]+32)+72)+40);
C[5]=+S[5]*(LambA*(LambA*(DkappA*(96*LambA+80*DkappA+144)+LambA*(16*LambA+320)+144)+
 DkappA*(48*(-DkappA-1)))+DkappA*(12*DkappA+24)+12);
C[4]=+S[5]*(LambA*(LambA*(16*(DkappA*(DkappA+1)+S[1]+1)+LambA*(32*DkappA+64))+DkappA*(
 16*(-DkappA-1)))+DkappA*(4*DkappA+8)+4);
C[3]=+32*S[8];
C[2]=+16*S[8];
C[1]=+4*S[2];
S[9]=V[8]*V[8];
S[9]*=S[9];
C[0]=+S[9];
}
