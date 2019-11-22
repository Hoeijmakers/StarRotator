#include <stdio.h>
#include <fcntl.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
long intdiv();
long mmin(),mmax();
float *vector();
double clight = 299791.0;
void convolv();
char *ggets(char *s);

main(int argc, char *argv[])
{
  int fd;
  int fo;
  double teff,logg,MH,vturb,start,dwave,wave;
  float Depth;
  char infile[80], ofile[80],tmp[80];
  float *depth;
  long i,k,n,npt,low,high,nspace,fs,num,nd;
  unsigned long left,qsize;
  double scale,vsini,space,sum1,sum2,a,z,Ds,n1,n2,ss;
  double l,d,u;
  double s2,st,dw;
  float *ys;
  char name[15],tmp1[10];
  int ni,nbytes;


  qsize = 80000;
  depth = vector(1,qsize);

  if(argc != 5) {
    printf("\nEnter name of input file > ");
    ggets(infile);
    printf("\nEnter name of output file > ");
    ggets(ofile);
    printf("\nEnter vsini (km/s) > ");
    ni = scanf("%lf",&vsini);
    printf("\nEnter limb darkening coefficient u > ");
    ni = scanf("%lf",&u);
    printf("\nSpacing in the output binary file will be the same\n");
    printf("as the input file\n");
    ggets(tmp);
  } else {
    strcpy(infile,argv[1]);
    strcpy(ofile,argv[2]);
    vsini = atof(argv[3]);
    u = atof(argv[4]);
  }


  if((fd = open(infile,O_RDWR)) == -1) {
    printf("\nError opening input file\n");
    exit(1);
  }
  

  nbytes = read(fd,&teff,sizeof(double));
  nbytes = read(fd,&logg,sizeof(double));
  nbytes = read(fd,&MH,sizeof(double));
  nbytes = read(fd,&vturb,sizeof(double));
  nbytes = read(fd,&start,sizeof(double));
  nbytes = read(fd,&dwave,sizeof(double));
  wave = start;
  k = 1;
  while(read(fd,&Depth,sizeof(float)) > 0) {
    depth[k] = Depth;
    k++;
  }
  close(fd);

  fo = open(ofile,O_CREAT|O_TRUNC|O_RDWR|S_IREAD|S_IWRITE,0666);
  nbytes = write(fo,&teff,sizeof(double));
  nbytes = write(fo,&logg,sizeof(double));
  nbytes = write(fo,&MH,sizeof(double));
  nbytes = write(fo,&vturb,sizeof(double));
  nbytes = write(fo,&start,sizeof(double));
  nbytes = write(fo,&dwave,sizeof(double));


  num = k-1;

  ys = vector(1,k-1);  


  s2 = (start+num*dwave)*vsini/(dwave*clight);
  nd = s2 + 5.5;

  convolv(depth,ys,num,nd,start,dwave,vsini,u);


  wave = start;
  for(i=1;i<=num;i++) nbytes = write(fo,&ys[i],sizeof(float));
  close(fo);

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
