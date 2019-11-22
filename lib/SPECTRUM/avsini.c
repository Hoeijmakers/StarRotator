#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
long intdiv();
long mmin(),mmax();
float *vector();
double clight = 299791.0;
void convolv();
char *ggets(char *s);

main(int argc, char *argv[])
{
  FILE *in, *out;
  double teff,logg,MH,vturb,start,dwave,wave;
  float Depth,wav;
  char infile[80], outfile[80],tmp[80];
  float *depth;
  long i,k,n,npt,low,high,nspace,fs,num,nd;
  unsigned long left,qsize;
  double scale,vsini,space,sum1,sum2,a,z,Ds,n1,n2,ss;
  double l,d,u;
  double s2,st,dw;
  float *ys;
  char name[15],tmp1[10];
  int ni;


  qsize = 200000;
  depth = vector(1,qsize);

  if(argc != 6) {
    printf("\nEnter name of input file > ");
    ggets(infile);
    printf("\nEnter name of output file > ");
    ggets(outfile);
    printf("\nEnter vsini (km/s) > ");
    ni = scanf("%lf",&vsini);
    printf("\nEnter limb darkening coefficient u > ");
    ni = scanf("%lf",&u);
    printf("\nEnter spacing in the input file; output spacing will be\n");
    printf("the same > ");
    ni = scanf("%lf",&dwave);
    ggets(tmp);
  } else {
    strcpy(infile,argv[1]);
    strcpy(outfile,argv[2]);
    vsini = atof(argv[3]);
    u = atof(argv[4]);
    dwave = atof(argv[5]);
  }


  if((in = fopen(infile,"r")) == NULL) {
    printf("\nError opening input file\n");
    exit(1);
  }
  
  out = fopen(outfile,"w");

  k = 1;
  while(fscanf(in,"%f %f",&wav,&depth[k]) != EOF) {
    if(k == 1) wave = wav;
    k++;
  }


  num = k-1;

  ys = vector(1,num);  
  start = wave;

  s2 = (start+num*dwave)*vsini/(dwave*clight);
  nd = s2 + 5.5;

  convolv(depth,ys,num,nd,start,dwave,vsini,u);


  for(i=1;i<=num;i++) {
    fprintf(out,"%9.3f %e\n",wave,ys[i]);
    wave += dwave;
  }
  fclose(out);

}

void convolv(y,ys,num,nd,st,dw,vsini,u)
float *y,*ys;
double st,dw,vsini,u;
long num,nd;
{
  double beta,gam,w,s,t,dlc,c1,c2,dv,r1,r2,f,v;
  long i,n1,n2,n;
  char tmp[10];

  beta = (1.0-u)/(1.0 - 0.333333*u);
  gam = u/(1.0 - 0.333333*u);  

/* End Effect */

  n1 = nd + 1;
  for(i=1;i<=nd;i++) ys[i] = y[i];
  n2 = num - nd -1;
  for(i=n2;i<=num;i++) ys[i] = y[i];
  if(vsini < 0.5) {
    for(i=1;i<=num;i++) ys[i] = y[i];
    return;
  }

/* convolve with rotation profile */

   w = st + (n1-1)*dw;
   for(n=n1;n<=n2;n++) {
     w = w+dw;
     s = 0.0;
     t = 0.0;
     dlc = w*vsini/clight;
     c1 = 0.63661977*beta/dlc;
     c2 = 0.5*gam/dlc;
     dv = dw/dlc;

     for(i=-nd;i<=nd;i++) {
       v = i*dv;
       r2 = 1.0 - v*v;
       if(r2 > 0.0) {
         f = c1*sqrt(r2) + c2*r2;
         t = t+f;
         s = s + f*y[n+i];
       }
     }
     ys[n] = s/t;
   }
   return;
}
