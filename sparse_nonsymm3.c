#include<stdio.h>
#include<stdlib.h>
int main(int argc,char *argv[]){
  int i, j, k, l, m, w;
  FILE *fp;
  fp = fopen(argv[2],"w");
  m=atoi(argv[1]);
  double a;
  fprintf(fp,"%d\n",m);
  fprintf(fp,"%d\n",m);
  fprintf(fp,"%d\n",2*m-1);
  for(i=0;i<m;i++){
      for(j=0;j<2;j++){
        if(i==m-1 && j==1) continue;
        if (j==0) a=1;
        else a=-1;
        fprintf(fp,"%d %d %30.20f\n",i,i+j,a);
    }
  }
  fclose(fp);
  return 0;
}
