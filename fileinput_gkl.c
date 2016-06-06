#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<sys/time.h>
double gettimeofday_sec(){
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return tv.tv_sec+(double)tv.tv_usec*1e-6;
}
#define TWO (2.0)
#define ONE (1.0)
#define ZERO (0.0)

int main(int argc, char *argv[]){

  char *matdescra="G  F  ";
  int izero=0,ione=1;
  double dzero=ZERO, done=ONE;

  int i, ii, j, lwork1, lwork2, m, max_bn, max_kk, max_kn, n, w, L ,K, accuracy;
  double alpha, scl, tmp, tmp2, t1, t2, sum0, sum1;
  char mode;
  L=atoi(argv[1]);
  accuracy = atoi(argv[3]);
  mode = argv[2][0];
  if(mode=='d') printf("dense mode だよ\n");
  else if(mode=='s') printf(",,,,,sparse mode だよ\n");
  else {printf("error: d[dense] か s[sparse] を指定して下さい。"); return 0;}

  K=L*2;
  scanf("%d",&m);
  scanf("%d",&n);
  scanf("%d",&w);

  printf("each row has %d elements.\n",w/m);
  int *indxA, *pntrbA, *pntreA;
  double *A;

  pntrbA=(int *)malloc( sizeof(int)*m );
  pntreA=(int *)malloc( sizeof(int)*m );

  if( pntrbA==NULL || pntreA==NULL){
    printf("Out of memory.\n");
    return 0;
  }

  indxA=(int *)malloc(sizeof(double)*w);
  if(indxA==NULL){
    printf("Out of memory.\n");
    return 0;
  }

  A=(double *)malloc(sizeof(double)*w);
  if(A==NULL){
    printf("Out of memory.\n");
    return 0;
  }
  
  if(mode=='s'){
    ii=0;
    w=0;
    pntrbA[0]=w+1;
    sum0=ZERO;
    while(scanf("%d %d %lf",&i,&j,&alpha) !=EOF){
      if(i!=ii){
        pntrbA[ii+1]=w+1;
        pntreA[ii]=pntrbA[ii+1];
        ii++;
      }
      indxA[w]=j+1;
      sum0=sum0+alpha*alpha;
      A[w]=alpha;
      w++;
    }
    pntreA[ii]=w+1;
  }
  else if(mode=='d'){
    w=0;
    while(scanf("%lf",&alpha) !=EOF){
      A[w]=alpha;
      w++;
    }
  }
  resgkl_main_(&mode,&accuracy,&m,&n,&L,&K,matdescra,indxA,pntrbA,pntreA,A);

  return 0;
}
