#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
long intdiv();
long mmin(),mmax();
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
  double l,d;
  char name[15];
  int ni;

  qsize = 1000000;
  depth = (float *) calloc(qsize,(unsigned long)sizeof(float));

  if(argc != 6) {
    printf("\nEnter name of input file > ");
    ggets(infile);
    printf("\nEnter name of output file > ");
    ggets(ofile);
    printf("\nEnter spacing in Angstroms of the input spectrum >");
    ni = scanf("%lf",&dwave);
    printf("\nEnter output resolution in Angstroms > ");
    ni = scanf("%lf",&res);
    printf("\nEnter spacing in Angstroms of spectrum in output file > ");
    ni = scanf("%lf",&space);
  } else {
    strcpy(infile,argv[1]);
    strcpy(ofile,argv[2]);
    dwave = atof(argv[3]);
    res = atof(argv[4]);
    space = atof(argv[5]);
  }

  if((in = fopen(infile,"r")) == NULL) {
    printf("\nError openint input file\n");
    exit(1);
  }

  fp = fopen(ofile,"w");

  k = 0;
  while(fscanf(in,"%f %f",&wav,&depth[k]) != EOF) {
    if(k == 0) wave = wav;
    k++;
  }


  npt = k - 1;

  res /= 2.0;
  n = intdiv(res,dwave);
  nspace = intdiv(space,dwave);
  a = -log10(0.50)/(n*n);

  i = 0;
  while(i <= npt) {
    low = (long) mmax(0,i-3*n);
    high = (long) mmin(npt,i+3*n);
    sum1 = sum2 = 0.0;
    for(k=low;k<=high;k++) {
      z = pow(10.0,-a*(k-i)*(k-i));
      sum1 += depth[k]*z;
      sum2 += z;
    }
    Ds = sum1/sum2;
    fprintf(fp,"%8.3f  %g\n",wave,Ds);
    wave += space;
    i += nspace;
  }
  free(depth);
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

long mmin(m,n)
long m,n;
{
  if(m <= n) return(m);
  else return(n);
}

long mmax(m,n)
long m,n;
{
  if(m >= n) return(m);
  else return(n);
}

