#include <math.h>
#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "spectrum.h"
double abund();
double dmax();
double dmin();
void qround();
int maxcharge();

void linelst(wave,list,nlist,atom,teff,logg,qf,reset)
double wave,teff,logg;
linelist *list;
atominfo *atom;
int *nlist,reset;
FILE *qf;
{
  int i,l,n;
  extern double inc;
  extern int flagr;
  extern int flagI;
  double lambda,code,El,Eu,loggf,ab,dampfac,SA,alp,sig;
  double gammar,gammas,gammaw,gam; 
  /* Gamma stark per electron number */
  /* Gamma van der waals per neutral hydrogen number */
  char tr[5];
  char tmp[10],buffer[100],src[20];
  static double lastwave = 0.0;
  double radmax = 10.0;
  double wave1,gfA,rad,lograd,neff,chi;
  double charg;
  double err;
  double hfsfac = 1.0;
  static int flag = 0;
  static int qfeof = 0;
  double fac = 1.0;
  int iso = 0;
  int flagai = 0;

  if(flag == 0 || reset == 1) lastwave = wave - radmax;

  if(wave < 3646.0) radmax = 20.0;
  else if(wave >= 8000.0) radmax = 20.0;
  else radmax = 10.0;
  wave1 = wave + radmax;

  if(qfeof == 1) return;
  if(lastwave > wave1) return;
  while(lastwave <= wave1) {
    if(fgets(buffer,100,qf) == NULL) {
      qfeof = 1;
      break;
    }
    /* printf("%s",buffer); */
    if(buffer[0] == '#') continue;
    lambda = atof(strtok(buffer," "));
    hfsfac = 1.0;
    code = atof(strtok(NULL," "));
    /* a provision for hyperfine splitting components */
    if(code < 0.0) {
      hfsfac = 3.0;
      code = fabs(code);
    }
    iso = 0;
    if(flagI == 1) iso = atoi(strtok(NULL," "));
    El = atof(strtok(NULL," "));
    Eu = atof(strtok(NULL," "));
    loggf = atof(strtok(NULL," "));
    dampfac = atof(strtok(NULL," "));
    strcpy(tr,strtok(NULL," "));
    if(strcmp(tr,"AI") == 0) flagai = 1;
    else flagai = 0;
    /* If transition type is AO, read in alp and sig parameters */
    if(strcmp(tr,"AO") == 0) {
      SA = atof(strtok(NULL," "));
      sig = floor(SA);
      alp = SA - floor(SA);
    } else alp = sig = 0.0;
    /* If transition type is GA, read in individual gammas, where
       Gamma stark is per electron number, Gamma van der Waals per
       neutral hydrogen number.  Logarithms of these Gammas should
       appear in the linelist */
    /* If the transition type is AI (autoionizing), gammar = Gamma_shore, 
       gammas = Ashore and gammaw = Bshore.  However, if Ashore, the asymmetry
       parameter, is negative, -log(|Ashore|), which will be a positive
       quantity, should be entered into the linelist file.  If Ashore is 
       positive, the entry is log(Ashore), which is negative.  The
       code below translates the various possibilities for GA and AI
       transition types */  
    if(strcmp(tr,"GA") == 0 || strcmp(tr,"AI") == 0) {
      gammar = pow(10.0,atof(strtok(NULL," ")));
      gam = atof(strtok(NULL," "));
      if(strcmp(tr,"AI") == 0) {
	if(gam >= 0.0) gammas = -pow(10.0,-gam);
	else gammas = pow(10.0,gam);
      } else gammas = pow(10.0,gam);
      gammaw = pow(10.0,atof(strtok(NULL," ")));
    } else gammar = gammas = gammaw = 0.0;
    err = atoi(strtok(NULL," "));
    if(err == 99) printf("You may need to use the i switch\n");
    if(lambda < lastwave) continue;
    lastwave = lambda;
    if(flagr == 1) qround(&lambda,inc);
/* Skip over molecules if effective temperature is too high */
    if(teff >= 8500.0 && code > 101.0) continue;
/* If atom is not supported by SPECTRUM, skip over the line */
    if(abund(atom,(int)floor(code)) < 0.0) continue;
/* If an unsupported ion of a supported atom is in the list, skip */
    if(maxcharge(atom,code) == -1) continue;
    list[*nlist].wave = lambda;
    list[*nlist].code = code;
    list[*nlist].iso = iso;
    list[*nlist].El = El;
    El = El*1.23981e-04;
    list[*nlist].Eu = Eu;
    list[*nlist].loggf = loggf;
    list[*nlist].fac = dampfac;
    strcpy(list[*nlist].T,tr);
    list[*nlist].alp = alp;
    list[*nlist].sig = sig;
    list[*nlist].gammar = gammar;
    list[*nlist].gammas = gammas;
    list[*nlist].gammaw = gammaw;
    if(flagai == 1) list[*nlist].ai = 1;
    else list[*nlist].ai = 0;
    if(flagai == 1) rad = radmax;
    else {
      ab = abund(atom,(int)floor(list[*nlist].code));
      gfA = loggf + log10(ab);
      charg = list[*nlist].code - floor(list[*nlist].code);
      if(list[*nlist].code - floor(list[*nlist].code) >= 0.01) lograd = 3.96 +
	   (0.75861 + 0.0258773*gfA)*gfA -0.242223*El;
      else lograd = 4.46 + (0.741125 + 0.0215996*gfA)*gfA +
	  ((0.0639498*El - 0.311767)*El + 0.0891479)*El;
      if(lograd > log10(radmax)) rad = radmax;
      else rad = pow(10.0,lograd);
      /* Make a rough temperature adjustment in the radius */
      if(teff > 5700.0) fac = 1.0;
      else fac = 1.0+(5700.0-teff)/400.0;
      rad *= fac;
      if(logg > 4.0) fac = 1.0+ logg-4.0;
      rad *= fac;
      /* Make adjustment for high excitation lines like Mg II 4481 */
      if(code < 100.0 && Eu > 7.0e+04) rad *= 6.0;
      /* Make adjustment for neutral lines near ionization limit */
      if(code - floor(code) < 0.05 && code < 100.0) {
        for(l=0;l<NATOM;l++) {
           if((int)floor(code) == atom[l].code) {
             i = l;
             break;
           }
        } 
        chi = atom[i].I1;
        if(chi > Eu*1.23981e-04) {
          neff = sqrt(13.585/(chi-Eu*1.23981e-04));
          if(neff > 10.0) rad *= pow(neff/10.0,2.0);
        }
      } 
      /* Make adjustment for CNO atomic lines */
      if(code >= 6.0 && code < 9.0) rad *= 4.0;
      rad = (double)dmax(rad,0.20);
      /* if(code > 100.0) rad = (double)dmin(rad,0.5); */
      /* Account for increasing doppler widths */
      if(lambda > 5000.0) rad *= lambda/5000.0;
      /* Make adjustment for minimum in H- opacity */
      if(lambda > 8000.0 && lambda <= 13000.0) rad *= 2.7;
      if(lambda > 13000.0 && lambda < 19000.0) rad *= 3.8;
      if(lambda >= 19000.0 && lambda < 50000.0) rad *= 2.7;
      /* Make adjustment for non-neutral species */
      if(code < 100 && charg > 0.0) 
        rad *= 3.5;
      if(code < 100 && charg > 0.1)
        rad *= 1.5;
      /* Increase computation radius for low-excitation lines */
      if(El < 0.124) rad *= 3.0;
      /* Increase computation radius for hyperfine components */
      rad *= hfsfac;
      /* Limit computation radius to radmax */
      rad = (double)dmin(rad,radmax);
    }
    list[*nlist].wavel = lambda - rad;
    list[*nlist].waveh = lambda + rad;
    list[*nlist].flag = 0;

    (*nlist)++;
    if(*nlist >= NLIST-2) break;
  }

  if(flag == 0) {
    flag = 1;
    return;
  }

  if(*nlist <= 1) return;
  while(list[0].wave <= wave) {
    for(n=0;n<(*nlist)-1;n++) {
       list[n].wave = list[n+1].wave;
       list[n].code = list[n+1].code;
       list[n].iso = list[n+1].iso;
       list[n].loggf = list[n+1].loggf;
       list[n].El = list[n+1].El;
       list[n].Eu = list[n+1].Eu;
       list[n].wavel = list[n+1].wavel;
       list[n].waveh = list[n+1].waveh;
       list[n].fac = list[n+1].fac;
       list[n].flag = list[n+1].flag;
       list[n].alp = list[n+1].alp;
       list[n].sig = list[n+1].sig;
       list[n].gammar = list[n+1].gammar;
       list[n].gammas = list[n+1].gammas;
       list[n].gammaw = list[n+1].gammaw;
       list[n].ai = list[n+1].ai;
       strcpy(list[n].T,list[n+1].T);
    }
    (*nlist)--;
    if(*nlist <= 1) break;
  }

}

