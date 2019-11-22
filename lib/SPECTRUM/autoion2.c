/* This routine is experimental only.  Don't believe anything in it! */

#include <stdio.h>
#include "spectrum.h"
#include <math.h>
struct autoionize {
  double code;
  double center;
  double rad;
  double Gam;
  double El;
  double Eh;
  double loggf;
};
double pfunction();
double fano(double eta);

double autoion(int ntau,atmosphere *model,double wave,pfunc *V)
{
  double kap,kappa = 0;
  double T,kT,U,eta;
  double Gcm;
  double e0,e;
  double El,Eh;
  double gf;
  double fac = 80.0;
  int N = 3;
  int i;
  static struct autoionize autol[3] = {
    {20.0,6318.109,15.0,1.0,35730.0,51554.0,0.699},
    {20.0,6343.308,15.0,1.1,35819.0,51579.0,0.845},
    {20.0,6361.786,15.0,1.5,35897.0,51611.0,0.954}};

  kT = model->kT[ntau];
  kappa = 0.0;
  for(i=0;i<N;i++) {
    if(wave < autol[i].center - autol[i].rad || 
       wave > autol[i].center + autol[i].rad) continue;

    Gcm = 1.0e+08*autol[i].Gam/(autol[i].center*autol[i].center);

    El = autol[i].El*1.23981e-04;
    Eh = autol[i].Eh*1.23981e-04;
    gf = pow(10.0,autol[i].loggf);

    e0 = 1.0e+08/autol[i].center;

    U = pfunction(V,autol[i].code,ntau);

    e = 1.0e+08/wave;
    eta = (e-e0)/Gcm;

    kap = pow((autol[i].center/wave),2.0)*(2.90e-06/autol[i].Gam)*
          fac*gf*fano(eta);
    kap *= 2.65386e-10*model->NCaI[ntau]*exp(-El/kT)*(1.0-exp((El-Eh)/kT))/U;
    kappa += kap;
  }
  return(kappa);
}

double fano(double eta)
{
  return(1.0/((eta*eta + 1.0)));
}
