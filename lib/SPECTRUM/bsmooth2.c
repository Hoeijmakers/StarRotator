
#include <stdio.h>
#include <fcntl.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
long intdiv();
long mmin(),mmax();
char *ggets(char *s);

main(int argc, char *argv[])
{
  int fd;
  FILE *fp;
  double teff,logg,MH,vturb,start,dwave,wave;
  float Depth;
  char infile[80], ofile[80],tmp[80];
  float *depth;
  long i,k,n,npt,low,high,nspace,fs;
  unsigned long left,qsize;
  double scale,res,space,sum1,sum2,a,z,Ds,n1,n2,ss;
  double l,d;
  char name[15];
  int ni;

  qsize = 1000000;
  depth = (float *) calloc(qsize,(unsigned long)sizeof(float));

  if(argc != 5) {
    printf("\nEnter name of input file > ");
    ggets(infile);
    printf("\nEnter name of output file > ");
    ggets(ofile);
    printf("\nEnter output resolution in Angstroms > ");
    ni = scanf("%lf",&res);
    printf("\nEnter spacing in Angstroms of spectrum in output file > ");
    ni = scanf("%lf",&space);
  } else {
    strcpy(infile,argv[1]);
    strcpy(ofile,argv[2]);
    res = atof(argv[3]);
    space = atof(argv[4]);
  }

  if((fd = open(infile,O_RDWR)) == -1) {
    printf("\nError opening input file\n");
    exit(1);
  }
  fp = fopen(ofile,"w");

  read(fd,&teff,sizeof(double));
  read(fd,&logg,sizeof(double));
  read(fd,&MH,sizeof(double));
  read(fd,&vturb,sizeof(double));
  read(fd,&start,sizeof(double));
  read(fd,&dwave,sizeof(double));
  wave = start;

  k = 0;
  while(read(fd,&Depth,sizeof(float)) > 0) {
    depth[k] = Depth;
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
    fprintf(fp,"%8.3f  %f\n",wave,1.0-Ds);
    wave += space;
    i += nspace;
  }
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

