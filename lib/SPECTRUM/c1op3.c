#include <math.h>
#include "spectrum.h"
#include <errno.h>
double partfn(double code, double T, double Ne);

double c1op_av(double nu,atmosphere *model,int j)
{

  double hckt,A,H,waveno,x,u;

  u = partfn(6.0,model->T[j],model->Ne[j]);
  waveno = 3.335641e-11*nu;
  hckt = 1.438832/model->T[j];

  H = 1.0e-30;

  /* 5921.789 A */
  if(waveno >= 16886.790) {
    x = 0.0;
    A = x*1.0*exp(-73975.91*hckt)*model->stim[j];
    H += A;
  }

  /* 5478.858 A */
  if(waveno >= 18251.980) {
    x = 0.0;
    A = x*5.0*exp(-72610.72*hckt)*model->stim[j];
    H += A;
  }

  /* 5131.416 A */
  if(waveno >= 19487.800) {
    x = 0.0;
    A = x*9.0*exp(-71374.90*hckt)*model->stim[j];
    H += A;
  }

  /* 4970.488 A */
  if(waveno >= 20118.750) {
    x = 0.0;
    A = x*3.0*exp(-70743.95*hckt)*model->stim[j];
    H += A;
  }

  /* 4730.212 A */
  if(waveno >= 21140.700) {
    x = 0.0;
    A = x*15.0*exp(-69722.00*hckt)*model->stim[j];
    H += A;
  }

  /* 4544.139 A */
  if(waveno >= 22006.370) {
    x = 2.1e-18*pow(22006.370/waveno,1.5);
    A = x*3.0*exp(-68856.33*hckt)*model->stim[j];
    H += A;
  }

  /* 3462.498 A */
  if(waveno >= 28880.880) {
    x = 1.54e-18*pow(28880.880/waveno,1.2);
    A = x*3.0*exp(-61981.82*hckt)*model->stim[j];
    H += A;
  }

  /* 3279.796 A */
  if(waveno >= 30489.700) {
    x = 0.2e-18*pow(30489.700/waveno,1.2);
    A = x*9.0*exp(-60373.00*hckt)*model->stim[j];
    H += A;
  }

  /* 1706.448 A */
  if(waveno >= 58601.270) {
    x = 0.0;
    A = x*9.0*exp(-75254.93*hckt)*model->stim[j];
    H += A;
  }

  /* 1445.663 A */
  if(waveno >= 69172.400) {
    x = 33.6e-18*pow(69172.400/waveno,1.5) -
	24.0e-18*pow(69172.400/waveno,2.5);
    x /= 3.0;
    A = x*1.0*exp(-21648.02*hckt)*model->stim[j];
    H += A;
  }

  /* 1444.339 A  */
  if(waveno >= 69235.820) {
    A *= 2.0;
    H += A;
  }

  /* 1433.335 A */
  if(waveno >= 69767.350) {
    x = 16.0e-18*pow(69767.350/waveno,3.0);
    A = x*15.0*exp(-64088.85*hckt)*model->stim[j];
    H += A;
  }

  /* 1240.268 A */
  if(waveno >= 80627.760) {
    x = 28.7e-18*pow(80627.760/waveno,1.5) -
	18.4e-18*pow(80627.760/waveno,2.5);
    x /= 3.0;
    A = x*5.0*exp(-10192.66*hckt)*model->stim[j];
    H += A;
  }

  /* 1239.293 A */
  if(waveno >= 80691.180) {
    A *= 2.0;
    H += A;
  }

  /* 1101.601 A */
  if(waveno >= 90777.00) {
    x = 38.6e-18*pow(90777.00/waveno,2.0) -
	28.2e-18*pow(90777.00/waveno,3.0);
    x /= 3.0;
    A = x*5.0*exp(-43.42*hckt)*model->stim[j];
    H += A;
  }

  /* 1101.273 A */
  if(waveno >= 90804.000) {
    A = x*3.0*exp(-16.42*hckt)*model->stim[j];
    H += A;
  }

  /* 1101.074 A */
  if(waveno >= 90820.420) {
    A = x*1.0*1.0*model->stim[j];
    H += A;
  }

  /* 1100.832 A */
  if(waveno >= 90840.420) {
    x *= 2.0;
    A = x*5.0*exp(-43.42*hckt)*model->stim[j];
    H += A;
  }

  /* 1100.504 A */
  if(waveno >= 90867.420) {
    A = x*3.0*exp(-16.42*hckt)*model->stim[j];
    H += A;
  }

  /* 1100.306 A */
  if(waveno >= 90883.840) {
    A = x*1.0*1.0*model->stim[j];
    H += A;
  }

  /* 998.792 A */
  if(waveno >= 100121.000) {
    x = 1.0e-18*pow(100121.000/waveno,3.0);
    A = x*5.0*exp(-33735.20*hckt)*model->stim[j];
    H += A;
  }

  return(H/u);
}








