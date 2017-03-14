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
  fprintf(fp,"%d\n",2*m);
  for(i=0;i<m;i++){
      for(j=0;j<2;j++){
        if(i==m-1 && j==1) { fprintf(fp,"%d %d %30.20f\n",i,i-1,0.0);
        continue;}
        fprintf(fp,"%d %d %30.20f\n",i,i+j,1.0);
    }
  }
  fclose(fp);
  return 0;
}
