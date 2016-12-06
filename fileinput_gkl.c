#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<sys/time.h>

#define TWO (2.0)
#define ONE (1.0)
#define ZERO (0.0)

int main(int argc, char *argv[]){

  int s,u,i,each, ii, j, m, n,tmp, w, L ,K, accuracy;
  int lwork, matrixseed, initseed, method;
  int maxthreads,k0,k1,t;
  double alpha,n0,l0,l1;
  char mode,ls;
  int *IAP,*JA,*start_row;
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
  maxthreads=omp_get_max_threads();
  start_row=(int *)malloc((maxthreads+1)*sizeof(int));
  if(IAP==NULL || JA==NULL || A==NULL || work==NULL || start_row==NULL){
    printf("Out of memory.\n");
    return 0;
  }

  n0 = (double)w/maxthreads;

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


  for(t=0,i=0;i<maxthreads && t<m;i++){
    start_row[i]=t;

    k1=0;
    while(k1<n0 && t<m){
      k0=k1;
      k1=k0+IAP[t+1]-IAP[t];
      t=t+1;
    }

    if(start_row[i]+1==t){
      continue;
    }

    l0=n0-k0;
    l1=fabs(k1-n0);
    if(l0<l1){
      t=t-1;
    }
  }
  start_row[i]=m;

  if(i!=maxthreads){
    for(i=i+1;i<maxthreads+1;i++){
      start_row[i]=0;
    }
  }
  for(i = 0;i<maxthreads+1;i++){
  printf("start_row %d = %d\n",i,start_row[i]);
  }
  
  resgkl_main_(start_row,&initseed,&method,&mode,&ls,&accuracy,&m,&n,&L,&K,IAP,JA,A,work,&lwork);

  return 0;
}

