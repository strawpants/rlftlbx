/* This is a collection of functions which  returns the permutation vector of a string array side1 which orders itself to side2 */
/* When parameters are present in side1 but not in side2 those will be put at the back or front of the permutation */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/* Define a new type holding a string address and 2 integers */

struct Acell{
  int a; /* Holds a tertiary sorting parameter */
  char* name;/*  Pointer to a string address */
};

/* Global parameter which hold the sizes of each of the sides */
int gn1,gn2;
int gst,gnd; /* Start and end for strng comparison */
int stlen; /* parameter string length */
size_t starti;
int ncom=0;
/* Function to return the FORTRAN index (starting with 1) based on the reference length */
size_t fortind(char *refer){
  return (size_t) (refer-starti)/stlen+1;
}

/* Function to return the FORTRAN index (starting with 1) based on a given reference length */
size_t fortind2(char *refer, size_t start){
  return (size_t) (refer-start)/stlen+1;
}

/* Comparison function which sorts alphabetically as the first criteria and storage wise as the second criteria
/* Returns a negative value if a1 should be in front of a2 and positive value if a1 should be behind a2 */
/* int comp_Alpha(struct Acell* a1, struct Acell* a2){ */
int comp_Alpha(const void* a, const void* b){
  int cmp;
/* Define temporay strings which hold copies of the valid part of the string with a null character appended */
  char str1[gnd-gst+2];
  char str2[gnd-gst+2];
  struct Acell a1=*(struct Acell*)a;
  struct Acell a2=*(struct Acell*)b;


  memcpy (str1, &a1.name[gst-1], (size_t)gnd-gst+1);
  memcpy (str2, &a2.name[gst-1], (size_t)gnd-gst+1);
  /* Add null characters at the end*/
  str1[gnd-gst+1]='\0';
  str2[gnd-gst+1]='\0';
  cmp=strcmp(str1,str2);

  if(cmp==0){ /* Then non-unique parameter is found */
    return (a1.a-a2.a); /* compare 2nd criteria */
  }else{
    return cmp;
  }
}

/* Second comnparison function (compare second criteria only) */
int comp_unique(const void* a, const void* b){
  struct Acell a1=*(struct Acell*)a;
  struct Acell a2=*(struct Acell*)b;
  return a1.a-a2.a;
}

/* Function to check whether the next entry is unique */
/* Returns 0 if unique and 1 if not  */
int check_unique(struct Acell * a1, struct Acell * a2){
  int cmp;
/* Define temporay strings which hold copies of the valid part of the string with a null character appended */
  char str1[gnd-gst+2];
  char str2[gnd-gst+2];
  memcpy (str1, &a1->name[gst-1], (size_t)gnd-gst+1);
  memcpy (str2, &a2->name[gst-1], (size_t)gnd-gst+1);
  /* Add null characters at the end (fortran doesn't do this)*/
  str1[gnd-gst+1]='\0';
  str2[gnd-gst+1]='\0';
  cmp=strcmp(str1,str2);
  if(cmp == 0){
    a1->a=a2->a-gn2;
    return 1;
  }else{
    return 0;
  }
}

/* Function to check whether the next entry is unique */
/* Returns 0 if unique and 1 if not  */
/* Check_unique2 also modifies the sorting parameter of the second side */
int check_unique2(struct Acell * a1, struct Acell * a2){
  int cmp;
/* Define temporay strings which hold copies of the valid part of the string with a null character appended */
  char str1[gnd-gst+2];
  char str2[gnd-gst+2];
  memcpy (str1, &a1->name[gst-1], (size_t)gnd-gst+1);
  memcpy (str2, &a2->name[gst-1], (size_t)gnd-gst+1);
  /* Add null characters at the end (fortran doesn't do this)*/
  str1[gnd-gst+1]='\0';
  str2[gnd-gst+1]='\0';
  cmp=strcmp(str1,str2);
  if(cmp == 0){
    a1->a=a2->a-2*gn2;
    a2->a=a2->a-gn2;
    return 1;
  }else{
    return 0;
  }
}


/* This function is to be called from fortran with a two contigous character arrays */
int get_permvecc_(int* n1, int* n2, int* slen,int* st,int* nd, int * back,int pvec[], char side1[][*slen], char side2[][*slen],int flen1,int flen2){
  int i,j;
/*   struct Acell dat[(*n1+*n2)]; /\* Create an array of Acell structures *\/ */
  struct Acell *dat = (struct Acell *)calloc(*n1+*n2,sizeof(struct Acell));
  /* Copy the values in the global sz1,sz2 parameters */
  gn1=*n1;
  gn2=*n2;
  gst=*st;
  gnd=*nd;
  stlen=*slen;

  starti=(size_t)&side1[0][0];

/* Fill up dat structure array */

     for(i=0 ;i < gn1 ;i++){
       dat[i].name=&side1[i][0]; /* Copy pointer value */
       dat[i].a=i-gn1-gn2;
    }
     for(i=0 ;i < gn2 ;i++){
       dat[i+gn1].name=&side2[i][0]; /* Copy pointer value */
       dat[i+gn1].a=i;
    }

     /* Sort array Alphabetically*/

     qsort(dat,(size_t)(gn1+gn2),sizeof(struct Acell),comp_Alpha); 


/* Now check for non-unique parameters */
     ncom=0;
     for(i=0;i<gn1+gn2-1;i++){
       if(dat[i].a < -gn2){
	 ncom+=check_unique(&dat[i],&dat[i+1]);
       }
     }
/* Apply second sorting */
     qsort(dat,(size_t)(gn1+gn2),sizeof(struct Acell),comp_unique); 
     
     
     /* Construct permutation vector */
     if(*back){ /* Put the unique parameters at the back */
       for(i=0;i<ncom;i++){
	 pvec[i]=fortind(dat[i+gn1-ncom].name);
       }
       for(i=ncom;i<gn1;i++){
	 pvec[i]=fortind(dat[i-ncom].name);
       }

     }else{
       for(i=0;i<gn1;i++){
	 pvec[i]=fortind(dat[i].name);
       }
     }
     free(dat);/*  Free up memory */
     return ncom;
}

/* Tis is the second versionb of the permve rouitne also allowing the second permutation vector to be retrieved */
/* The back parmeter is removed (one can exchange the pvec1 and pvec2 to achieve the same result) */
/*To be called from fortran with a two contigous character arrays */
int get_permvecc2_(int* n1, int* n2, int* slen,int* st,int* nd,int * back, int pvec1[],int pvec2[], char side1[][*slen], char side2[][*slen],int flen1,int flen2){
  int i,j;
/*   struct Acell dat[(*n1+*n2)]; /\* Create an array of Acell structures *\/ */
  struct Acell *dat = (struct Acell *)calloc(*n1+*n2,sizeof(struct Acell));
  /* Copy the values in the global sz1,sz2 parameters */
  gn1=*n1;
  gn2=*n2;
  gst=*st;
  gnd=*nd;
  stlen=*slen;

 size_t start1=(size_t)&side1[0][0];;
 size_t start2=(size_t)&side2[0][0];;
/* Fill up dat structure array */

     for(i=0 ;i < gn1 ;i++){
       dat[i].name=&side1[i][0]; /* Copy pointer value */
       dat[i].a=i-gn1-gn2;
    }
     for(i=0 ;i < gn2 ;i++){
       dat[i+gn1].name=&side2[i][0]; /* Copy pointer value */
       dat[i+gn1].a=i+gn2;
    }

     /* Sort array Alphabetically ( and on memory order as second criteria)*/

     qsort(dat,(size_t)(gn1+gn2),sizeof(struct Acell),comp_Alpha); 

/* Now check for non-unique parameters */
     ncom=0;
     for(i=0;i<gn1+gn2-1;i++){
       if(dat[i].a < -gn2){
	 ncom+=check_unique2(&dat[i],&dat[i+1]);
       }
     }
     
/* Apply second sorting */
     qsort(dat,(size_t)(gn1+gn2),sizeof(struct Acell),comp_unique); 

/*      for(i=0;i<gn1+gn2;i++){ */
/*        printf("after %d %d\n",i,dat[i].a); */
/*      } */
     


     /* Construct permutation vector */
     if(*back){ /* Put the unique parameters at the back for permutation1 */
       for(i=0;i<ncom;i++){
	 pvec1[i]=fortind2(dat[i+gn1-ncom].name,start1);
       }
       for(i=ncom;i<gn1;i++){
	 pvec1[i]=fortind2(dat[i-ncom].name,start1);
       }

     }else{
       for(i=0;i<gn1;i++){
	 pvec1[i]=fortind2(dat[i].name,start1);
       }
     }

     for(i=0;i<gn2;i++){
       pvec2[i]=fortind2(dat[i+gn1].name,start2);
     }
     free(dat);/*  Free up memory */
     return ncom;
}


