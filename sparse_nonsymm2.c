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
  fprintf(fp,"%d\n",3*m-2);
  for(i=0;i<m;i++){
      for(j=0;j<3;j++){
        if(i==0 && j==0) continue;
        if(i==m-1 && j==2) continue;
        if (i==0 && j==1) a=1;
        else if(j==1) a=2;
        else a=-1;
        fprintf(fp,"%d %d %30.20f\n",i,i-1+j,a);
    }
  }
  fclose(fp);
  return 0;
}
