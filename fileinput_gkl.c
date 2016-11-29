#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<sys/time.h>

#define TWO (2.0)
#define ONE (1.0)
#define ZERO (0.0)

int main(int argc, char *argv[]){

  int s,u,i,each, ii, j, m, n,tmp, w, L ,K, accuracy, lwork, matrixseed, initseed, method;
  double alpha;
  char mode,ls;
  int *IAP,*JA;
  double *A,*work;

  FILE *fp;
  method = atoi(argv[1]);//1=qr1,2=qr2,3=oqds1,4=oqds2
  L=atoi(argv[2]);
  accuracy = atoi(argv[5]);
  mode = argv[3][0];
  ls = argv[4][0];
  matrixseed = atoi(argv[6]);
  initseed = atoi(argv[7]);
  fp = fopen(argv[8],"r");
  if(mode=='d') printf("runnning: dense mode\n");
  else if(mode=='s') printf("runnning: sparse mode\n");
  else {printf("error: d[dense] か s[sparse] を指定して下さい。"); return 0;}

  if(ls=='s') printf("runnning: smallest singularvalues ... \n");
  if(ls=='l') printf("runnning: largest singularvalues ... \n");
  if(method==1) printf("runnning: 1...QR(biside singular vectors) ... \n");
  if(method==2) printf("runnning: 2...QR(oneside singular vectors) ... \n");
  if(method==3) printf("runnning: 3...OQDS(biside singular vectors) ... \n");
  if(method==4) printf("runnning: 4...OQDS(oneside singular vectors) ... \n");

  K=L*2;
  fscanf(fp,"%d",&m);
  fscanf(fp,"%d",&n);
  fscanf(fp,"%d",&w);

  printf("each row has %d elements.\n",w/m);

  srand(matrixseed);
 
  IAP=(int *)malloc(sizeof(int)*(m+1));
  JA=(int *)malloc(sizeof(int)*w);
  A=(double *)malloc(sizeof(double)*w);

  lwork=5*K;
  if (m>lwork) lwork=m;
  if (n>lwork) lwork=n;
  if (K*K>lwork) lwork=K*K;

  work=(double *)malloc(sizeof(double)*lwork);

  if(IAP==NULL || JA==NULL || A==NULL || work==NULL){
    printf("Out of memory.\n");
    return 0;
  }

  each = w / m;
  if(mode=='s'){
    ii=0;
    IAP[0]=1;
    for(tmp=0;tmp<w;tmp++){
    //while(fscanf("%d %d %lf",&i,&j,&alpha) !=EOF){
      i = tmp / each;
      j = rand() % n;
      alpha = rand()/(double)RAND_MAX;
      //printf("%d %d %30.20f\n",i,j,alpha);
      if(i!=ii){
        IAP[ii+1]=tmp+1;
        ii++;
      }
      JA[tmp]=j+1;
      A[tmp]=alpha;
      //tmp++;
    }
    IAP[ii+1]=tmp+1;
  }
  else if(mode=='d'){
    tmp=0;
    while(fscanf(fp,"%lf",&alpha) !=EOF){
      A[tmp]=alpha;
      tmp++;
    }
  }
  fclose(fp);
  resgkl_main_(&initseed,&method,&mode,&ls,&accuracy,&m,&n,&L,&K,IAP,JA,A,work,&lwork);

  return 0;
}

