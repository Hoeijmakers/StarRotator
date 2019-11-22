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
void hotDensity(atmosphere *model,atominfo *atom,double ah,double ahe,int flagw);
void veryhotDensity(atmosphere *model,atominfo *atom,double ah,double ahe,int flagw);
void tauref();
void tauwave();
double flux();
int lline();
void pop();
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
void inisotope();
void isorelabun();
memo reset;
void setreset();
int flagk = 0;
int flagmgh = 0;
int flagp = 0;
int flagP = 0;
int flagt = 0;
int flagI = 0;
int flagu = 0;
int flagC = 0;
int flagcd = 0;
int flagO = 0;
int flags = 0;
int flagSq = 0;
int mghla = 0;
int mghlb = 0;
int NI = 0;
char *ggets(char *s);
/* variables for isotopes */
double ra1H,ra2H,ra12C,ra13C,ra14N,ra15N,ra16O,ra17O,ra18O;
double ra24Mg,ra25Mg,ra26Mg,ra28Si,ra29Si,ra30Si,ra40Ca,ra42Ca;
double ra43Ca,ra44Ca,ra46Ca,ra48Ca,ra46Ti,ra47Ti,ra48Ti,ra49Ti;
double ra50Ti;
FILE *opout;
char buf2[100];


int main(int argc, char *argv[])
{
  int i,j,k,flag;
  int flagv = 0;
  int nflag = 0;
  long n;
  int code;
  int flagw = 1;
  double ah,ahe,waveref,wave,Flux,Depth,w,w0,vturb,vt,ew,original_abund;
  double nabund,abund0,abund1,abundmid,wmid,Atot,AH;
  double start, end;
  double eqmin = 0.0;
  atmosphere *model;
  atominfo *atom;
  linelist *list;
  linedata *line;
  isodata *isotope;
  pfunc *V;
  population *POP;
  char *file,*ofile,name[60],oname[60],lines[60],*flines,c,atmdat[60];
  char tmp[10],isofile[80],select[80];
  FILE *qf,*fp,*sel;
  int nq = 0;
  int ni;

  setreset(0);
  strcpy(atmdat,"stdatom.dat");
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
  if((isotope = (isodata *) calloc(500,sizeof(isodata))) == NULL)
    nrerror("Allocation of memory for isotope failed");

  printf("\nLINES v2.75b  (C) Richard O. Gray May 13, 2008\n\n");

  if(argc > 1) {
    if(strcspn(argv[1],"v") != strlen(argv[1])) flagv = 1;
    if(strcspn(argv[1],"i") != strlen(argv[1])) flagI = 1;
    if(strcspn(argv[1],"t") != strlen(argv[1])) flagt = 1;
    if(strcspn(argv[1],"C") != strlen(argv[1])) flagC = 1;
    if(strcspn(argv[1],"s") != strlen(argv[1])) flags = 1;
  }
  flagcd = flagC;
  flagC = 0;

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
  if(flagI == 1) {
    printf("\nEnter name of isotope data file (default = isotope.iso) > ");
    ggets(isofile);
    if(strcmp(isofile,"") == 0) strcpy(isofile,"isotope.iso");
    inisotope(isofile,isotope);
    isorelabun(isotope);
  }
  printf("\nEnter starting, ending wavelength > ");
  ni = scanf("%lf,%lf",&start,&end);
  ggets(tmp);
  printf("\nEnter name of output file > ");
  ofile = ggets(oname);
  fp = fopen(ofile,"w");

  printf("\nEnter name of atom data file Default = stdatom.dat > ");
  ggets(atmdat);
  if(strcmp(atmdat,"") == 0) strcpy(atmdat,"stdatom.dat");

  if(flags == 1) {
    printf("\nEnter name of output file for selected lines > ");
    ggets(select);
    if((sel = fopen(select,"w")) == NULL) {
      printf("\nCannot open output file for selected lines\n");
      exit(1);
    }
  }

  printf("\nEnter vturb (km/s)> ");
  ni = scanf("%lf",&vt);
  vt *= 1.0e+05;
  for(i=0;i<Ntau;i++) model->mtv[i] = vt;

  printf("\nMinimum Equivalent Width for output (mA) > ");
  ni = scanf("%lf",&eqmin);
  
  ah = 0.911;
  ahe = 0.089;
  inatom(atmdat,atom,model->MH,&ah,&ahe);


  pfinit(V,atom,model,flagw);
  printf("\nCalculating Number Densities\n");
  if(model->teff <= 8500) Density(model,atom,ah,ahe,flagw);
  else if(model->teff >= 25000.0) veryhotDensity(model,atom,ah,ahe,flagw); 
  else hotDensity(model,atom,ah,ahe,flagw);
  popinit(POP,atom,model,V,flagw);

  waveref = 5000.0;
  printf("Calculating Reference Opacities\n");
  tauref(model,waveref);
  printf("Entering Main Loop\n");

  while((nq = lline(&wave,line,atom,qf,start,end,isotope)) >= 1) {
    if(wave < start) continue;
    if(wave > end) break;
    if(nq == 2) continue;
    if(nq == 3) continue;
    flag = 0;
    /*    vturb = vt*1.0e+05; */
    Depth = 1.0;
    ew = eqwidth(model,line,wave,V,POP,&Depth);
    if(ew < eqmin && flagSq == 0) continue;
    fprintf(fp,"%8.3f  %5.1f  %8.3f  %6.4f\n",line[0].wave,line[0].code,
	       ew,1.0-Depth);
    fclose(fp);
    if(flags == 1) fprintf(sel,"%s",buf2);
    if(flagv == 1) printf("%8.3f  %5.1f  %8.3f  %6.4f\n",
                              line[0].wave,line[0].code,ew,1.0-Depth);
    fp = fopen(ofile,"a");
  }
  fclose(fp);
  if(flags == 1) fclose(sel);
  return(0);

}

double eqwidth(model,line,wave,V,POP,Depth)
atmosphere *model;
linedata *line;
double wave;
pfunc *V;
population *POP;
double *Depth;
{
  double dwave = 0.005;
  double w1,w2,w,d0,Flux;
  int n = 0;
  int flagq;
  double q;

  flagq = 0;
  tauwave(model,wave);
  Flux = flux(model,wave);
  pop(line,0,model,V,POP);
  broad(model,line,0,line[0].sig,line[0].alp,line[0].fac);
  capnu(line,0,model);

  w = w1 = w2 = 0.0;
  while(1) {
    eqtaukap(wave,model,line);
    if(flagq == 0 && flagcd == 1) flagC = 1;
    w2 = depth(model,wave,Flux);
    if(flagq == 0) {
      *Depth = w2;
      flagq = 1;
    }
    if(flagq == 1 && flagcd == 1) flagC = 0;
    if(n == 0) d0 = w2;
    if(n > 0) w += 0.5*(w1 + w2)*dwave;
    w1 = w2;
/*  if(w2 < 0.10*d0) dwave = 0.01;
    if(w2 < 0.02*d0) dwave = 0.05;  */
    if(w2 < 0.00001) return(2000*w);
    wave += dwave;
    n++;
  }
}





