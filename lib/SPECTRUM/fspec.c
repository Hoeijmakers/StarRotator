#include <stdio.h>
#include "spectrum.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
int Ntau = 72;
void inmodel();
void Density();
void tauref();
void tauwave();
void tauflx();
double flux();
void fluxflx();
double intensity();
void inlin();
void pop();
void broad();
void capnu();
void taukap();
void inatom();
void linelst();
void pfinit();
void popinit();
double depth();
double depthflx();
double depthmu();
float **cmatrix();
void nrerror();
float **bkap;
double interval();
memo reset;
void setreset();
double dmax(),dmin();
int imax(),imin();
double inc;
int flagr = 0;
int flagc = 0;
int flagk = 0;
int flagg = 0;
float *velgrad;
double mu = 1.0;

void fspec(char *name, char *lines, char *atmdat, char *oname, double vturb, 
double mu, double start, double end, double dwave, char *flag)
{
  int i,k;
  int nline = 0;
  int nlist = 0;
  int code;
  int flagf = 0;
  int flagm = 0;
  int flagi = 1;
  int flagb = 0;
  int flagw = 1;
  int flaga = 0;
  int flagd = 0;
  int flagM = 0;
  int flagG = 0;
  double ah,ahe,waveref,start,end,wave,dwave,dw,Flux,Depth,wstart,wend,vturb;
  double nabund,DW=20.0,Intensity,flxstart,flxend,q,y;
  float qd;
  unsigned len_written;
  atmosphere *model;
  atominfo *atom;
  linelist *list;
  linedata *line;
  linedata *strgln;
  pfunc *V;
  population *POP;
  Helium *He;
  char *file,name[80],oname[80],lines[80],c,tmp[20],atmdat[80];
  char *vg,vgrad[80];
  char *ofile;
  FILE *qf,*fp;
  int fd;
  char junk[10];


  setreset(0);
  strcpy(atmdat,"atom.dat");

  if(argc > 1) {
    if(strcspn(flag,"f") != strlen(flag)) flagf = 1;
    if(strcspn(flag,"m") != strlen(flag)) flagm = 1;
    if(strcspn(flag,"M") != strlen(flag)) flagM = 1;
    if(strcspn(flag,"b") != strlen(flag)) flagb = 1;
    if(strcspn(flag,"n") != strlen(flag)) flagw = 0;
    if(strcspn(flag,"a") != strlen(flag)) flaga = 1;
    if(strcspn(flag,"d") != strlen(flag)) flagd = 1;
    if(strcspn(flag,"r") != strlen(flag)) flagr = 1;
    if(strcspn(flag,"c") != strlen(flag)) flagc = 1;
    if(strcspn(flag,"k") != strlen(flag)) flagk = 1;
    if(strcspn(flag,"g") != strlen(flag)) flagg = 1;
    if(strcspn(flag,"G") != strlen(flag)) flagG = 1;
  }
  flaga = 1;
  flagw = 0;
  flagG = flagg = 0;
  if(flagm == 1 || flagM == 1) flagf = 0;
  if(flagm == 1 || flagf == 1 || flagM == 1) flagi = 0;
  if(flagf == 1 && flagr == 0)
    printf("\nPlease consider using switch r with switch f\n");

  if((model = (atmosphere *) calloc(1,sizeof(atmosphere))) == NULL)
    nrerror("Allocation of memory for atmosphere failed");
  if((V = (pfunc *) calloc(1,sizeof(pfunc))) == NULL)
    nrerror("Allocation of memory for pfunc failed");
  if((atom = (atominfo *) calloc(NATOM,sizeof(atominfo))) == NULL)
    nrerror("Allocation of memory for atom failed");
  if((POP = (population *) calloc(NATOM-NMOL,sizeof(population))) == NULL)
    nrerror("Allocation of memory for population failed");
  if((line = (linedata *)
  calloc((unsigned long)NLINE,(unsigned long)sizeof(linedata))) == NULL)
    nrerror("Allocation of memory for line failed");

  if((strgln = (linedata *)
    calloc(NSTRG,(unsigned long)sizeof(linedata))) == NULL)
    nrerror("Allocation of memory for caii failed");
  if((list = (linelist *) calloc(NLIST,sizeof(linelist))) == NULL)
    nrerror("Allocation of memory for list failed");
  if((He = (Helium *) calloc(NHE,sizeof(Helium))) == NULL)
    nrerror("Allocation of memory for He failed");
  if((velgrad = (float *) calloc(NTAU,sizeof(float))) == NULL)
    nrerror("Allocation of memory for velgrad failed");
  bkap = cmatrix(0,3,0,NTAU);

  
  file = *name;
  inmodel(model,file,flagw);
  if((qf = fopen(lines,"r")) == NULL) {
   printf("Cannot find line data file\n");
   exit(1);
  }
  ofile = *oname;
  vturb *= 1.0e+05;
  inc = dwave;

  if(flagb == 0) fp = fopen(ofile,"w");
  if(flagb == 1) {
    fd = open(oname,O_CREAT|O_TRUNC|O_RDWR|S_IREAD|S_IWRITE);
    if(fd == -1) {
     perror("Error:");
     exit(1);
    }
    write(fd,&model->teff,sizeof(double));
    write(fd,&model->logg,sizeof(double));
    write(fd,&model->MH,sizeof(double));
    q = vturb/1.0e+05;
    write(fd,&q,sizeof(double));
    write(fd,&start,sizeof(double));
    write(fd,&dwave,sizeof(double));
  }
  /* Normal abundances for hydrogen and helium; inatom may change these */
  ah = 0.911;
  ahe = 0.089;
  inatom(atmdat,atom,model->MH,&ah,&ahe);

  pfinit(V,atom,model,flagw);
  nline = 0;
  nlist = 0;
  Density(model,atom,ah,ahe,flagw);
  popinit(POP,atom,model,V,flagw);

  waveref = 5000.0;
  tauref(model,waveref);
  wstart = start;
  DW = interval(start,dwave,flagf);
  wend = start + DW;
  wend = (double)dmin(wend,end);

  while(wstart < end - 0.0001) {
     wave = (wstart + wend)/2.0;
     if(flagi == 1 || flagm == 1 || flagM == 1) tauwave(model,wave);
     if(flagf == 1) tauflx(model,wstart,wend);
     if(flagi == 1) Flux = flux(model,wave);
     if(flagf == 1) fluxflx(model,wstart,wend,&flxstart,&flxend);
     if(flagm == 1 || flagM == 1) Intensity = intensity(model,wave,mu);

     wave = wstart;

     while(wend-wave > dwave/2.0 + 0.0001) {
	linelst(wave,list,&nlist,atom,model->teff,qf,0);
	inlin(wave,line,&nline,list,nlist,atom);
	for(i=0;i<nline;i++) {
	  if(line[i].flag == 0) {
	  pop(line,i,model,V,POP);
	  broad(model,line,i,vturb,line[i].fac);
	  capnu(line,i,model);
	  }
	}
	taukap(wave,model,atom,line,nline,strgln,vturb,V,He,POP);
	if(flagi == 1) Depth = depth(model,wave,Flux);
	if(flagf == 1) Depth = depthflx(model,wave,wstart,wend);
	if(flagm == 1 || flagM == 1) Depth = depthmu(model,wave,Intensity,mu);

	if(flagb == 0 && (flagi == 1 || flagm == 1))
	       fprintf(fp,"%9.3f   %f\n",wave,1.0 - Depth);
        if(flagb == 0 && flagM == 1)
               fprintf(fp,"%9.3f   %e\n",wave,Intensity*(1.0 - Depth));
	if(flagb == 1 && (flagi == 1 || flagm == 1)) {
	       qd = Depth;
	       write(fd,&qd,sizeof(float));
	}
        if(flagb == 1 && flagM == 1) {
               qd = Intensity*(1.0-Depth);
               write(fd,&qd,sizeof(float));
        }
	if(flagf == 1) {
	  Flux = flxstart + (flxend-flxstart)*(wave-wstart)/(wend-wstart);
          if(Depth > Flux) Depth = Flux;
	  if(flagb == 0)
		   fprintf(fp,"%9.3f   %e\n",wave,Flux*(1.0-Depth/Flux));
	  if(flagb == 1) {
		   qd = Flux*(1.0-Depth/Flux);
		   write(fd,&qd,sizeof(float));
	  }
	}
	wave += dwave;
     }

     wstart = wave;
     DW = interval(wave,dwave,flagf);
     DW = (double)dmax(dwave,DW);
     wend = wstart + DW;
     wend = (double) dmin(wend,end);

  }
  fclose(qf);
  if(flagb == 0) fclose(fp);
  if(flagb == 1) close(fd);
}
