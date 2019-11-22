#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

main(int argc, char *argv[])
{
  char buffer[150],*tmp,buf2[150];
  char infile[80],outfile[80];
  double teff,logg,mh;
  double TEFF,LOGG,MH;
  FILE *fp,*out;
  int i,N;

  if(argc != 5) {
    printf("\nPlease use the following format:\n");
    printf("\nselectmod supermodel.dat output.mod teff logg\n\n");
    exit(1);
  }

  printf("\nWarning:  if this program results in a segmentation fault,\n");
  printf("it means that in the supermodel you must translate <CR> to\n");
  printf("<CR><LF>\n"); 

  strcpy(infile,argv[1]);
  strcpy(outfile,argv[2]);
  TEFF = atof(argv[3]);
  LOGG = atof(argv[4]);

  if((fp = fopen(infile,"r")) == NULL) {
    printf("\nCannot find input supermodel\n");
    exit(1);
  }
  out = fopen(outfile,"w"); 


  while(1) {
    do {
      if(fgets(buffer,120,fp) == NULL) {
	printf("File access error\n");
	exit(1);
      }
    } while(strstr(buffer,"TEFF") == NULL);
    strcpy(buf2,strstr(buffer,"TEFF"));
    tmp = strtok(buf2," ");
    teff = atof(strtok(NULL," "));
    if(teff != TEFF) continue;
    strcpy(buf2,strstr(buffer,"GRAVITY"));
    tmp = strtok(buf2," ");
    logg = atof(strtok(NULL," "));
    if(logg != LOGG) continue;
    fprintf(out,"%s",buffer);
    for(i=0;i<96;i++) {
      if(fgets(buffer,120,fp) == NULL) {
	printf("File access error\n");
	exit(1);
      }
      fprintf(out,"%s",buffer);
    }
    fclose(fp);
    fclose(out);
    break;
  }
}

