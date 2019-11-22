#include <math.h>
#include <errno.h>
#include <stdio.h>
#include "spectrum.h"
double pfunction();
int approx();
void popinit();

/* H2 added Dec 13, 2007 */

void pop(line,N,model,V,POP)
linedata *line;
atmosphere *model;
pfunc *V;
population *POP;
int N;
{
  double Q4,ratio,cneutral,c1ion,c2ion,c3ion,mu,D0,Ua,Ub,Na,Nb,theta,psip;
  int i,charge;
  static int l = 2;
  double k = 8.617084e-05;
  double gffac = 1.0;
  extern int Ntau;

  if(approx(line[N].code,6.0,0.001) == 1) {
    for(i=0;i<Ntau;i++) line[N].xnum[i] = model->NCI[i]*
			exp(-line[N].El/(k*model->T[i]))/pfunction(V,6.0,i);
    return;
  }

  if(approx(line[N].code,7.0,0.001) == 1) {
    for(i=0;i<Ntau;i++) line[N].xnum[i] = model->NNI[i]*
			exp(-line[N].El/(k*model->T[i]))/pfunction(V,7.0,i);
    return;
  }

  if(approx(line[N].code,8.0,0.001) == 1) {
    for(i=0;i<Ntau;i++) line[N].xnum[i] = model->NOI[i]*
			exp(-line[N].El/(k*model->T[i]))/pfunction(V,8.0,i);
    return;
  }


  if(line[N].code < 100.0) {
    cneutral = floor(line[N].code);
    c1ion = cneutral + .1;
    c2ion = cneutral + .2;
    c3ion = cneutral + .3;
    if(approx(line[N].code,cneutral,0.001) == 1) charge = 0;
    else if(approx(line[N].code,c1ion,0.001) == 1) charge = 1;
    else if(approx(line[N].code,c2ion,0.001) == 1) charge = 2;
    else if(approx(line[N].code,c3ion,0.001) == 1) charge = 3;
    else {
	  printf("Could not determine charge in pop\n");
	  printf("Setting charge = 0\n");
	  charge = 0;
    }
    if(approx(POP[l].atom,line[N].code,0.3) != 1) {
      for(i=2;i<NATOM-NMOL;i++) {
	if(approx(POP[i].atom,line[N].code,0.3) == 1) {
	  l = i;
	  break;
	}
      }
    }

    for(i=0;i<Ntau;i++) {
      Q4 = POP[l].Q4[i];
      if(charge == 0) line[N].xnum[i] = line[N].abund*model->NA[i]*
	 exp(-1.1605e+4*line[N].El/model->T[i])/(Q4*pfunction(V,cneutral,i));
      if(charge == 1) line[N].xnum[i] = line[N].abund*model->NA[i]*
	 POP[l].R1[i]*exp(-1.1605e+4*line[N].El/model->T[i])/
	 (Q4*pfunction(V,c1ion,i));
      if(charge == 2) line[N].xnum[i] = line[N].abund*model->NA[i]*
	 POP[l].R1[i]*POP[l].R2[i]*exp(-1.1605e+4*line[N].El/model->T[i])/
	 (Q4*pfunction(V,c2ion,i));
      if(charge == 3) line[N].xnum[i] = line[N].abund*model->NA[i]*
	 POP[l].R1[i]*POP[l].R2[i]*POP[l].R3[i]*
	 exp(-1.1605e+4*line[N].El/model->T[i])/(Q4*pfunction(V,c3ion,i));
    }
    return;
  } else {
    mu = line[N].chi2;
    D0 = line[N].chi1;
    gffac = line[N].chi3;

    for(i=0;i<Ntau;i++) {
      if(approx(line[N].code,101.0,0.001) == 1) {
        Ua = pfunction(V,1.0,i);
        Ub = pfunction(V,1.0,i);
        Na = model->NHI[i];
        Nb = model->NHI[i];
      }
      if(approx(line[N].code,106.0,0.001) == 1) {
	Ua = pfunction(V,1.0,i);
	Ub = pfunction(V,6.0,i);
	Na = model->NHI[i];
	Nb = model->NCI[i];
      }
      if(approx(line[N].code,107.0,0.001) == 1) {
	Ua = pfunction(V,1.0,i);
	Ub = pfunction(V,7.0,i);
	Na = model->NHI[i];
	Nb = model->NNI[i];
      }
      if(approx(line[N].code,108.0,0.001) == 1) {
	Ua = pfunction(V,1.0,i);
	Ub = pfunction(V,8.0,i);
	Na = model->NHI[i];
	Nb = model->NOI[i];
      }
      if(approx(line[N].code,112.0,0.001) == 1) {
	Ua = pfunction(V,1.0,i);
	Ub = pfunction(V,12.0,i);
	Na = model->NHI[i];
	Nb = model->NMgI[i];
      }
      if(approx(line[N].code,114.0,0.001) == 1) {
	Ua = pfunction(V,1.0,i);
	Ub = pfunction(V,14.0,i);
	Na = model->NHI[i];
	Nb = model->NSiI[i];
      }
      if(approx(line[N].code,120.0,0.001) == 1) {
	Ua = pfunction(V,1.0,i);
	Ub = pfunction(V,20.0,i);
	Na = model->NHI[i];
	Nb = model->NCaI[i];
      }
      if(approx(line[N].code,606.0,0.001) == 1) {
	Ua = pfunction(V,6.0,i);
	Ub = pfunction(V,6.0,i);
	Na = model->NCI[i];
	Nb = model->NCI[i];
      }
      if(approx(line[N].code,607.0,0.001) == 1) {
	Ua = pfunction(V,6.0,i);
	Ub = pfunction(V,7.0,i);
	Na = model->NCI[i];
	Nb = model->NNI[i];
      }
      if(approx(line[N].code,608.0,0.001) == 1) {
	Ua = pfunction(V,6.0,i);
	Ub = pfunction(V,8.0,i);
	Na = model->NCI[i];
	Nb = model->NOI[i];
      }
      if(approx(line[N].code,814.0,0.001) == 1) {
	Ua = pfunction(V,8.0,i);
	Ub = pfunction(V,14.0,i);
	Na = model->NOI[i];
	Nb = model->NSiI[i];
      }
      if(approx(line[N].code,822.0,0.001) == 1) {
	Ua = pfunction(V,8.0,i);
	Ub = pfunction(V,22.0,i);
	Na = model->NOI[i];
	Nb = model->NTiI[i];
      }

      theta = 5040.0/model->T[i];
      psip = 1.38054e-16*model->T[i]*pow(10.0,D0*theta-13.670)*
	     pow(theta,2.5)/(pow(mu,1.5)*Ua*Ub);

      line[N].xnum[i] = psip*gffac*Na*Nb*exp(-line[N].El/(k*model->T[i]));
    }
    return;
  }

}

void popinit(POP,atom,model,V,flagw)
population *POP;
atominfo *atom;
atmosphere *model;
pfunc *V;
int flagw;
{
  int i,j;
  double Q1,Q2,Q3,R4,T15;
  extern int Ntau;

  if(flagw == 1) {
    printf("\nCalculating Ionization ratios for all atoms at all levels\n");
    printf("Completed atomic number: ");
  }
  for(i=2;i<NATOM-NMOL;i++) {
    POP[i].atom = atom[i].code;
    for(j=0;j<Ntau;j++) {
      Q1 = Q2 = Q3 = 1.0;
      T15 = pow(model->T[j],1.5);
      if(atom[i].maxcharge == 3)  {
	 R4 = (4.830e+15/model->Ne[j])*(pfunction(V,POP[i].atom+0.4,j)/
	 pfunction(V,POP[i].atom+0.3,j))*T15*
	 exp(-1.16e+4*atom[i].I4/model->T[j]);
	 Q1 += R4;
      }
      if(atom[i].maxcharge >= 2)  {
	 POP[i].R3[j] = (4.830e+15/model->Ne[j])*(pfunction(V,POP[i].atom+0.3,j)/
	 pfunction(V,POP[i].atom+0.2,j))*T15*
	 exp(-1.16e+4*atom[i].I3/model->T[j]);
	 Q2 += POP[i].R3[j]*Q1;
      }
      if(atom[i].maxcharge >= 1)  {
	 POP[i].R2[j] = (4.830e+15/model->Ne[j])*(pfunction(V,POP[i].atom+0.2,j)/
	 pfunction(V,POP[i].atom+0.1,j))*T15*
	 exp(-1.16e+4*atom[i].I2/model->T[j]);
	 Q3 += POP[i].R2[j]*Q2;
	 POP[i].R1[j] = (4.830e+15/model->Ne[j])*(pfunction(V,POP[i].atom+0.1,j)/
	 pfunction(V,POP[i].atom,j))*T15*exp(-1.16e+4*atom[i].I1/model->T[j]);
	 POP[i].Q4[j] = 1.0 + POP[i].R1[j]*Q3;
      }
    }
    if(flagw == 1) printf("%2.0f  ",POP[i].atom);
  }
  if(flagw == 1) printf("\n");
}


