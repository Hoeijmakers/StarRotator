#include <stdio.h>
#include "spectrum.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fcntl.h>
#define FACTOR 1.6
#define JMAX 100
int Ntau = 72;
void inmodel();
void Density();
void tauref();
void tauwave();
double flux();
int bwline();
void pop();
void broad(atmosphere *model,linedata *line,int N,double Sig,
           double Alp,double fac);
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
int flagI = 0;
int flagp = 0;
int flagP = 0;
int flagt = 0;
int flagu = 0;
int flagO = 0;
int flagC = 0;
int flagCNO = 0;
char *ggets(char *s);
FILE *opout;

main(int argc, char *argv[])
{
  int i,j,k,flag;
  int code;
  int flagw = 1;
  double ah,ahe,waveref,wave,Flux,Depth,w,w0,vturb,vt,ew,original_abund;
  double nabund,abund0,abund1,abundmid,wmid,Atot,AH,vtl,vth,vts;
  atmosphere *model;
  atominfo *atom;
  linelist *list;
  linedata *line;
  pfunc *V;
  population *POP;
  char *file,*ofile,name[60],oname[60],lines[60],*flines,c,atmdat[60];
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
  if((line = (linedata *) calloc(1,sizeof(linedata))) == NULL)
    nrerror("Allocation of memory for line failed");

  printf("\nBLACKWELL v2.75 (C) Richard O. Gray 1992, 2008\n\n");


  printf("\nEnter name of stellar atmosphere data file > ");
  file = ggets(name);
  inmodel(model,file,flagw);
  printf("\nEnter name of line input file > ");
  flines = ggets(lines);
  if(strcmp(lines,"") == 0) strcpy(lines,"bw.lst");
  if((qf = fopen(lines,"r")) == NULL) {
   printf("Cannot find line data file\n");
   exit(1);
  }
  printf("\nEnter name of output file > ");
  ofile = ggets(oname);
  fp = fopen(ofile,"w");

  printf("\nEnter name of atom data file Default = stdatom.dat > ");
  ggets(atmdat);
  if(strcmp(atmdat,"") == 0) strcpy(atmdat,"stdatom.dat");

  printf("\nEnter vturb low, high, step > ");
  ni = scanf("%lf,%lf,%lf",&vtl,&vth,&vts);

  
  ah = 0.911;
  ahe = 0.089;
  inatom(atmdat,atom,model->MH,&ah,&ahe);


  pfinit(V,atom,model,flagw);
  printf("\nCalculating Number Densities\n");
  Density(model,atom,ah,ahe,flagw);
  popinit(POP,atom,model,V,flagw);

  waveref = 5000.0;
  printf("Calculating Reference Opacities\n");
  tauref(model,waveref);
  printf("Entering Main Loop\n");

  while(bwline(&wave,line,atom,&ew,qf) == 1) {
     if(line[0].code == 6.0 || line[0].code == 7.0 || line[0].code == 8.0) 
       flagCNO = 1;
     else flagCNO = 0;
     original_abund = line[0].abund;
     flag = 0;
     for(vt=vtl;vt<=vth;vt+=vts) {
       vturb = vt*1.0e+05;
       for(i=0;i<Ntau;i++) model->mtv[i] = vturb;
       line[0].abund = original_abund;
       while(eqwidth(model,line,atom,wave,V,POP) > ew) {
	 line[0].abund /= FACTOR;
       }
       abund0 = line[0].abund;
       line[0].abund = original_abund;
       while(eqwidth(model,line,atom,wave,V,POP) < ew)
	   line[0].abund *= FACTOR;
       abund1 = line[0].abund;
       for(j=1;j<=JMAX;j++) {
	 abundmid = line[0].abund = (abund0 + abund1)/2.0;
	 wmid = eqwidth(model,line,atom,wave,V,POP);
	 if(wmid >= ew) abund1 = abundmid;
	 if(wmid < ew) abund0 = abundmid;
	 if(fabs(log10(abund1/abund0)) <= 0.0002) break;
       }
       Atot = log10((abund1+abund0)/2.0);
       AH = Atot + 0.0405 + 12.0;
       printf("%8.3f  %5.1f  %4.2f  %7.3f  %7.3f\n",line[0].wave,line[0].code,
	       vt,Atot,AH);
       fprintf(fp,"%8.3f  %5.1f  %4.2f  %7.3f  %7.3f\n",line[0].wave,line[0].code,
	       vt,Atot,AH);
       fclose(fp);
       fp = fopen(ofile,"a");
     }
  }

}


double eqwidth(model,line,atom,wave,V,POP)
atmosphere *model;
linedata *line;
atominfo *atom;
double wave;
pfunc *V;
population *POP;
{
  double dwave = 0.005;
  double w1,w2,w,d0,Flux;
  int i,m;
  int n = 0;

  if(flagCNO == 1) {
    m = 0;
    if(line[0].code == 6.0) m = 5;
    if(line[0].code == 7.0) m = 6;
    if(line[0].code == 8.0) m = 7;
    if(m == 0) {
      printf("\nCNO error in Blackwel ... now exiting\n");
      exit(1);
    }
  }
  tauwave(model,wave);
  Flux = flux(model,wave);
  pop(line,0,model,V,POP);
  /* Note that in pop8, for codes 6.0, 7.0 & 8.0, line[0].xnum is hardwired,
     and so eqwidth returns the same equivalent width for CI, NI and OI lines
     without the correction below.  Correction added Nov 7, 2006 */
  if(flagCNO == 1) {
    for(i=0;i<Ntau;i++) line[0].xnum[i] *= line[0].abund/atom[m].abund;
  }
  broad(model,line,0,line[0].sig,line[0].alp,line[0].fac);
  capnu(line,0,model);

  w = w1 = w2 = 0.0;
  while(1) {
    eqtaukap(wave,model,line);
    w2 = depth(model,wave,Flux);
    if(n == 0) d0 = w2;
    if(n > 0) w += 0.5*(w1 + w2)*dwave;
    w1 = w2;
    /* if(w2 < 0.10*d0) dwave = 0.01;
       if(w2 < 0.02*d0) dwave = 0.05; */
    if(w2 < 0.00001) return(2000*w);
    wave += dwave;
    n++;
  }
}
