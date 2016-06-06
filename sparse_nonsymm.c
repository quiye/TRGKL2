#include<stdio.h>
#include<stdlib.h>
int main(int argc,char *argv[]){
  int i, j, k, l, m, n, w;
  FILE *fp;
  fp = fopen(argv[4],"w");
  m=atoi(argv[1]);
  n=atoi(argv[2]);
  w=atoi(argv[3]);
  srand(1);
  double a;
  fprintf(fp,"%d\n",m);
  fprintf(fp,"%d\n",n);
  fprintf(fp,"%d\n",w*m);
  for(i=0;i<m;i++){
    for(j=0;j<w;j++){
      k=rand() % n;
      a=rand()/(double)RAND_MAX;
      fprintf(fp,"%d %d %30.20f\n",i,k,a);
    }
  }
  fclose(fp);
  return 0;
}
