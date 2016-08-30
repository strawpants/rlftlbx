/* selection of functions which performs numerical sorting using the qsort from the standard C lib */
/* These functions are intended to be called from fortran */
/* 15 June 2011:: RR added the double sorting functions */

#include <stdio.h>
#include <stdlib.h>

int inverse;
struct Ptrcell{
  int* pnt; /* pointer holding memory addresses of the orinal memory (unchanged) */
};

struct Ptrdcell{
  double* pnt; /* pointer holding memory addresses of the orinal memory (unchanged) */
};



int compare_int(const void *a, const void *b) /* Compare based on numerical order and secondly on position in memory */
{
  struct Ptrcell a1=*(struct Ptrcell*)a;
  struct Ptrcell a2=*(struct Ptrcell*)b;
  int c;
  c=*a1.pnt-*a2.pnt;
  if(inverse)c=-c;
  if(c==0){ /* Apply second criteria (memory location) */
    return a1.pnt-a2.pnt;
  }else{
    return c;
  }
}

int compare_dbl(const void *a, const void *b) /* Compare based on numerical order and secondly on position in memory */
{
  struct Ptrdcell a1=*(struct Ptrdcell*)a;
  struct Ptrdcell a2=*(struct Ptrdcell*)b;
  int c;
  c=*a1.pnt-*a2.pnt;
  if(inverse)c=-c;
  if(c==0){ /* Apply second criteria (memory location) */
    return a1.pnt-a2.pnt;
  }else{
    return c;
  }
}


void isort_c_(int pvec[],int *n,int *inv,int a[]){
  int i,j;
  struct Ptrcell *dat = (struct Ptrcell *)calloc(*n,sizeof(struct Ptrcell));

  inverse=*inv; /* set inverse option */
/*   copy memory locations in structure */
  for(i=0;i<*n;i++){
    dat[i].pnt=&a[i];
  }
  
/* sort the set of pointers */
  qsort(dat,(size_t)*n,sizeof(struct Ptrcell),compare_int);

  /* Retrieve the associated permutation vector */
  for(i=0;i<*n;i++){
/*     fprintf(stderr,"%d %p %p %d\n",*dat[i].pnt,dat[i].pnt,&a[0],(int)(dat[i].pnt-&a[0]+1)); */
    pvec[i]=(dat[i].pnt-&a[0])+1;
  }
  free(dat);
}

void dsort_c_(int pvec[],int *n,int *inv, double a[]){
  int i,j;
  struct Ptrdcell *dat = (struct Ptrdcell *)calloc(*n,sizeof(struct Ptrdcell));

  inverse=*inv; /* set inverse option */
/*   copy memory locations in structure */
  for(i=0;i<*n;i++){
    dat[i].pnt=&a[i];
  }
  
/* sort the set of pointers */
  qsort(dat,(size_t)*n,sizeof(struct Ptrdcell),compare_dbl);

  /* Retrieve the associated permutation vector */
  for(i=0;i<*n;i++){
/*     fprintf(stderr,"%d %p %p %d\n",*dat[i].pnt,dat[i].pnt,&a[0],(int)(dat[i].pnt-&a[0]+1)); */
    pvec[i]=(dat[i].pnt-&a[0])+1;
  }
  free(dat);
}







