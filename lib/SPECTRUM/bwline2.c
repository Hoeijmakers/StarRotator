#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "spectrum.h"
#include <string.h>
int approx();

int bwline(wave,line,atom,ew,qf)
double *wave,*ew;
linedata *line;
atominfo *atom;
FILE *qf;
{
  int i,icode,j;
  double lambda,code,code1,El,Eu,loggf,df,eqw,sig,alp,SA;
  double gammar,gammas,gammaw; 
  char T[5];
  char tmp[10],buffer[120];

  if(fgets(buffer,100,qf) == NULL) return(0);
  if(buffer[0] == '#') return(3);
  *wave = lambda = atof(strtok(buffer," "));
  line[0].wave = lambda;
  line[0].code = code = atof(strtok(NULL," "));
  line[0].El = 1.23981e-04*atof(strtok(NULL," "));
  line[0].Eu = 1.23981e-04*atof(strtok(NULL," "));
  loggf = atof(strtok(NULL," "));
  line[0].gf = pow(10.0,loggf);
  line[0].fac = atof(strtok(NULL," "));
  strcpy(T,strtok(NULL," "));
  strcpy(line[0].T,T);
  /* If transition type is AO, read in alp and sig parameters */
  if(strcmp(T,"AO") == 0) {
    SA = atof(strtok(NULL," "));
    sig = floor(SA);
    alp = SA - floor(SA);
  } else alp = sig = 0.0;
  /* If transition type is GA, read in individual gammas, where
       Gamma stark is per electron number, Gamma van der Waals per
       neutral hydrogen number.  Logarithms of these Gammas should
       appear in the linelist */
  if(strcmp(T,"GA") == 0) {
    gammar = pow(10.0,atof(strtok(NULL," ")));
    gammas = pow(10.0,atof(strtok(NULL," ")));
    gammaw = pow(10.0,atof(strtok(NULL," ")));
  } else gammar = gammas = gammaw = 0.0;
  line[0].alp = alp;
  line[0].sig = sig;
  line[0].gammar = gammar;
  line[0].gammas = gammas;
  line[0].gammaw = gammaw;
  *ew = atof(strtok(NULL," "));
  if(approx(code,floor(code),0.001) == 1) icode = 0;
  else if(approx(code,floor(code)+0.1,0.001) == 1) icode = 1;
  i = 0;
  while((int)code != atom[i].code) i++;
  line[0].atomass = atom[i].amass;
  line[0].chi1 = atom[i].I1;
  line[0].chi2 = atom[i].I2;
  if(icode == 0) line[0].chi = atom[i].I1;
  else line[0].chi = atom[i].I2;
  line[0].abund = atom[i].abund;
  line[0].flag = 0;
  for(j=0;j<NTAU;j++) line[0].xnum[j] = line[0].a[j] =
		      line[0].dopp[j] = 0.0;
  return(1);
}
