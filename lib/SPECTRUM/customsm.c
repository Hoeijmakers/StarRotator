#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
long intdiv();
double linespread();
int mmin();
int mmax();
char *ggets(char *s);

main(int argc, char *argv[])
{
  FILE *in;
  FILE *fp;
  double teff,logg,MH,vturb,start,dwave,wave;
  float Depth;
  char infile[80], ofile[80],tmp[80];
  float *depth,wav;
  long i,k,n,npt,low,high,nspace,fs;
  unsigned long left,qsize;
  double scale,res,space,sum1,sum2,a,z,Ds,n1,n2,ss;
  double l,d,x,w,norm,M,q;
  char name[15];
  double c = 2.772588722;
  double pi = 3.141592653589793;
  int ni;

  qsize = 320000;
  depth = (float *) calloc(qsize,(unsigned long)sizeof(float));

  if(argc != 7) {
    printf("\nEnter name of input file > ");
    ggets(infile);
    printf("\nEnter name of output file > ");
    ggets(ofile);
    printf("\nEnter spacing in Angstroms of the input spectrum > ");
    ni = scanf("%lf",&dwave);
    printf("\nEnter output resolution in Angstroms > ");
    ni = scanf("%lf",&w);
    printf("\nEnter spacing in Angstroms of spectrum in output file > ");
    ni = scanf("%lf",&space);
    printf("\nEnter shape factor between 0 and 1 > ");
    ni = scanf("%lf",&M);
  } else {
    strcpy(infile,argv[1]);
    strcpy(ofile,argv[2]);
    dwave = atof(argv[3]);
    w = atof(argv[4]);
    space = atof(argv[5]);
    M = atof(argv[6]);
  }
   
  norm = 0.5*M*pi*w + (1.0-M)*w*sqrt(pi/c);

  if((in = fopen(infile,"r")) == NULL) {
    printf("\nError openint input file\n");
    exit(1);
  }

  fp = fopen(ofile,"w");

  k = 0;
  while(fscanf(in,"%f %f",&wav,&depth[k]) != EOF) {
    printf("%f\n",wav);
    if(k == 0) wave = wav;
    if(k == 0) start = wave;
    k++;
  }


  npt = k - 1;

  while(wave <= start+dwave*npt) {
    n = intdiv(wave - start,dwave);
    low = (long)mmax(0,n-intdiv(10.0,dwave));
    high = (long)mmin(npt,n+intdiv(10.0,dwave));
    sum1 = sum2 = 0.0;
    for(k=low;k<=high;k++) {
      x = (start + k*dwave) - wave;
      q = linespread(x,M,w);
      sum1 += q*depth[k]*dwave;
      sum2 += q*dwave;
    }
    Ds = sum1/sum2;
    fprintf(fp,"%8.3f  %g\n",wave,Ds);
    wave += space;
  }
  fclose(in);
  fclose(fp);
}


long intdiv(a,b)
double a,b;
{
  double n1,n2,ab;
  long n;

  ab = a/b;
  n1 = floor(ab);
  n2 = ceil(ab);
  if(fabs(ab-n1) < fabs(ab-n2)) n = n1;
  else n = n2;

  return(n);
}

double linespread(X,M,w)
double X,M,w;
{
  double arg,lor,ex,y;
  double c = 2.772588722;  

  arg = X/w;
  ex = exp(-c*arg*arg);
  lor = 1.0/(4.0*arg*arg + 1.0);
  y = (M*lor + (1-M)*ex);
  return(y);

}

int mmin(m,n)
int m,n;
{
  if(m <= n) return(m);
  else return(n);
}

int mmax(m,n)
int m,n;
{
  if(m >= n) return(m);
  else return(n);
}
