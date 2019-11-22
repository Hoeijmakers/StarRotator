#include <stdio.h>
#include "spectrum.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fcntl.h>
#define FACTOR 1.6
#define JMAX 40
int Ntau = 72;
void inmodel();
void Density();
void tauref();
void tauwave();
double flux();
int bwline();
void pop();
void broad();
void broad(atmosphere *model,linedata *line,int N,double sig,
           double alp,double fac);
void capnu();
void eqtaukap();
void inatom();
void pfinit();
void popinit();
double depth();
void nrerror();
double abund();
double eqwidth();
memo reset;
void setreset();
int flagk = 0;
int flagmgh = 0;
int mghla = 0;
int mghlb = 0;
int flagp = 0;
int flagP = 0;
int flagt = 0;
int flagu = 0;
int flagC = 0;
int flagO = 0;
char *ggets(char *s);
FILE *opout;

main(int argc, char *argv[])
{
  int i,j,k,flag;
  int code;
  int flagw = 1;
  double ah,ahe,waveref,wave,Flux,Depth,w,w0,vturb,vt,ew,original_gf;
  double nabund,gf0,gf1,gfmid,wmid,loggf;
  atmosphere *model;
  atominfo *atom;
  linelist *list;
  linedata *line;
  pfunc *V;
  population *POP;
  char *file,*ofile,name[60],oname[60],lines[60],*flines,c,atmdat[60];
  char tmp[10];
  FILE *qf,*fp;
  int ni;


  setreset(0);
  strcpy(atmdat,"stdatom.dat");

  if(argc > 1) {
     if(strcspn(argv[1],"t") != strlen(argv[1])) flagt = 1;
  }

  if((model = (atmosphere *) calloc(1,sizeof(atmosphere))) == NULL)
    nrerror("Allocation of memory for atmosphere failed");
  if((V = (pfunc *) calloc(1,sizeof(pfunc))) == NULL)
    nrerror("Allocation of memory for pfunc failed");
  if((atom = (atominfo *) calloc(NATOM,sizeof(atominfo))) == NULL)
    nrerror("Allocation of memory for atom failed");
  if((POP = (population *) calloc(NATOM-NMOL,sizeof(population))) == NULL)
    nrerror("Allocation of memory for population failed");
  if((line = (linedata *)
    calloc(1,(unsigned long)sizeof(linedata))) == NULL)
    nrerror("Allocation of memory for line failed");

  printf("\nGFADJUST, (C) Richard O. Gray 1992, 2001\n\n");

  printf("\nEnter name of stellar atmosphere data file > ");
  file = ggets(name);
  inmodel(model,file,flagw);
  printf("\nEnter name of input MOD file > ");
  flines = ggets(lines);
  if(strcmp(lines,"") == 0) strcpy(lines,"bw.lst");
  if((qf = fopen(lines,"r")) == NULL) {
   printf("Cannot find line data file\n");
   exit(1);
  }
  printf("\nEnter name of output ADJ file > ");
  ofile = ggets(oname);

  printf("\nEnter microturbulence (km/s) > ");
  ni = scanf("%lf",&vturb);
  vturb *= 1.0e+05;
  for(i=0;i<Ntau;i++) model->mtv[i] = vturb;
  ggets(tmp);

  printf("\nEnter name of atom data file Default = atom.dat > ");
  ggets(atmdat);
  if(strcmp(atmdat,"") == 0) strcpy(atmdat,"atom.dat");

  fp = fopen(ofile,"w");
  
  ah = 0.911;
  ahe = 0.089;
  inatom(atmdat,atom,model->MH,&ah,&ahe);
  /* printf("\nWould you like to change any of the elemental abundances y/n 
? > ");
  c = getch();
  printf("\n");
  while(c == 'y' || c == 'Y') {
     printf("\nEnter atomic number of element > ");
     ni = scanf("%d",&code);
     for(i=0;i<NATOM;i++) {
       if(atom[i].code == code) k = i;
     }
     printf("\nPlease enter new abundance for element %d.",atom[k].code);
     printf("\nCurrent log(N/Ntotal) = %f > ",log10(atom[k].abund));
     ni = scanf("%lf",&nabund);
     atom[k].abund = pow(10.0,nabund);
     printf("\nWould you like to change another elemental abundance y/n ? > ");
     c = getch();
     printf("\n");
  }
  */

  pfinit(V,atom,model,flagw);

  printf("Calculating Number Densities\n");
  Density(model,atom,ah,ahe,flagw);
  popinit(POP,atom,model,V,flagw);

  waveref = 5000.0;
  printf("Calculating Reference Opacities\n");
  tauref(model,waveref);
  printf("Entering Main Loop\n");

  while(bwline(&wave,line,atom,&ew,qf) == 1) {
     tauwave(model,wave);
     Flux = flux(model,wave);
     pop(line,0,model,V,POP);
     broad(model,line,0,line[0].sig,line[0].alp,line[0].fac);
     original_gf = line[0].gf;
     flag = 0;
     line[0].gf = original_gf;
       while(eqwidth(model,line,wave,Flux) > ew)
	   line[0].gf /= FACTOR;
       gf0 = line[0].gf;
       line[0].gf = original_gf;
       while(eqwidth(model,line,wave,Flux) < ew)
	   line[0].gf *= FACTOR;
       gf1 = line[0].gf;
       for(j=1;j<=JMAX;j++) {
	 gfmid = line[0].gf = (gf0 + gf1)/2.0;
	 wmid = eqwidth(model,line,wave,Flux);
	 if(wmid >= ew) gf1 = gfmid;
	 if(wmid < ew) gf0 = gfmid;
	 if(fabs(log10(gf1/gf0)) <= 0.01) break;
       }
       loggf = log10((gf1+gf0)/2);
       printf("%8.3f  %5.1f  %5.2f\n",line[0].wave,line[0].code,
	       loggf);
       fprintf(fp,"%8.3f  %5.1f  %5.2f\n",line[0].wave,line[0].code,
	       loggf);
       fclose(fp);
       fp = fopen(ofile,"a");
  }

}




double eqwidth(model,line,wave,Flux)
atmosphere *model;
linedata *line;
double wave,Flux;
{
  double dwave = 0.005;
  double w1,w2,w,d0;
  int n = 0;


  capnu(line,0,model);

  w = w1 = w2 = 0.0;
  while(1) {
    eqtaukap(wave,model,line);
    w2 = depth(model,wave,Flux);
    if(n == 0) d0 = w2;
    if(n > 0) w += 0.5*(w1 + w2)*dwave;
    w1 = w2;
    if(w2 < 0.00001) return(2000*w);
    wave += dwave;
    n++;
  }
}



