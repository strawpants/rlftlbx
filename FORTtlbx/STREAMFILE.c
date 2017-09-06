/*C methods which open stream files for unformatted I/O
These routines are constructed such that they can be called from Fortran through a fortran wrapper around C
The idea is that these routines provide a portable implementation of Stream access files which may be piped on linux*/

//Coded by Roelof Rietbroek 15-07-2008
//updated 11-08-2009: catch errors for fread and fwrite
//updated 25-11-2009: added new swapping routine and duplicate functions without underscore (allows for different linker standards)
//Updated 09-02-2011 added fseek functionality
#include <stdio.h>
#include <stdlib.h>
//Declare a list of 300 file pointers
const int maxf=300;
FILE *filelist[300]; //Maximum of 300 files to be opened at the same time
int tracker[300] = {0}; //tracker which keeps hold of the openstatus of the files (initialize to unopened)

//Function to open a stream for writing opr reading
int opencstream_(int *perm,char *filename, long int flen){
  int unit=-1;
  int i=0;
  int err=0;
/*   fprintf(stderr, "Length of string %li\n",flen); */
/*   fprintf(stderr, "inputted File name %s\n",filename); */

  //Search for an open spot
  while(unit <1){
    if(i > maxf){
      fprintf(stderr,"opencstream: ERROR too many files open, exiting\n");
      exit(1);
    }
    if(tracker[i] == 0){

      unit=i+1;
      tracker[i]=1;
      
    }else{ 
      i++;
    }
  }

  switch(*perm){
  case 0: // open for  read only
    filelist[unit-1]=fopen(filename,"rb");
    break;
  case 1: //open for write only 
    filelist[unit-1]=fopen(filename,"wb");
    break;
  case 3: //read from standard input
    filelist[unit-1]=stdin;
    break;
  case 4: //write to standard output
    filelist[unit-1]=stdout;
    break;
  }
  if(filelist[unit-1]==NULL){
    fprintf(stderr,"opencstream: ERROR opening file, exiting\n");
    exit(1);
  }
  
/*   /\*Set buffer size to 8192 *\/ */

/*   err=setvbuf(filelist[unit-1],NULL,_IOFBF,8192); */
/*   if(err != 0){ */
/*     fprintf(stderr,"opencstream: ERROR setting buffer, exiting\n"); */
/*     exit(1); */
/*   } */
  return unit;
  
}



//function to close a file
void cclose_(int *unit){
  int err;
  //  fprintf(stderr,"unit %i \n",*unit);
  err=fclose(filelist[*unit-1]);
    if(err != 0){
      fprintf(stderr,"ERROR closing file, exiting \n");
      exit(1);
    }

  tracker[*unit-1]=0;
  
}


//Function to write arbitrary data to a file (integer, doubles, real, characters, etc)
// The actual data is the last argument so that it also allowed to pass characters arrays from fortran
void cwrite_( int *unit,int *bytes, const void *input){
    if(fwrite(input,1,*bytes,filelist[*unit-1]) != *bytes){
          fprintf(stderr,"ERROR in cwrite: incomplete write\n");
	  exit(1);
    }
  
}

//function which reads in arbitrary data from a file

void cread_( int *unit,int *bytes, void *output){
  if(fread(output,1,*bytes,filelist[*unit-1]) != *bytes){
    fprintf(stderr,"ERROR in cread: incomplete read\n");
    exit(1);
  }

}

void creadmany_( int *unit, size_t *bytes, void *output){
    //fprintf(stderr,"Bytes requested %ld\n",*bytes);
  if(fread(output,1,*bytes,filelist[*unit-1]) != *bytes){
    fprintf(stderr,"ERROR in cread: incomplete read\n");
    exit(1);
  }

}

//Function to skip forward in a file. Will return an error if the file is not seekable
void cskip_(int *unit,int *bytes){
  if(fseek(filelist[*unit-1],*bytes,SEEK_CUR) != 0){
    fprintf(stderr,"ERROR in cskip: file not seekable? (pipe/fifo?)\n");
    exit(1);
  }
}

//debug function which prints the memory adress of a pointer to standard error
void printaddress_(void *input){
  fprintf(stderr,"%p\n",input);
}

/* void printchar_(char *input){ */
/*   fprintf(stderr,"%s\n",input[1]); */
/* } */

/* /\* Function to copy a pointer adress *\/ */
/* void copyaddress32_(int* out, int* in){ */
/* /\*   fprintf(stdout,"in,out, %d %d \n",*in,*out); *\/ */
/*   *out=*in; */
/* /\*   fprintf(stdout,"in,out, %d %d \n",*in,*out); *\/ */
/* } */

/* /\* Function to copy a pointer adress *\/ */
/* void copyaddress64_(size_t* in,size_t* out, long int flen){ */
/*     fprintf(stdout,"in,out, %p %p \n",in,out); */
/*     out=in; */
/*   fprintf(stdout,"in,out, %p %p \n",in,out); */
/* } */

void cswap_(int * bytes, int * sz, unsigned char* data){
  int half=*bytes/2;
  unsigned char dum;
  int shft;
  int i,j;
  for( i=0 ; i< *sz ; i++){
    shft=*bytes *i;
     for(j=0; j< half;j++){
       dum=data[shft+j];
       data[shft+j]=data[shft+*bytes-j-1];
       data[shft+*bytes-j-1]=dum;
            } 
  }

}



/* BELOW ARE THE SAME ROUTINES WITHOUT UNDERSCORES ADDED FOR COMPATIBILITY */
 
//Function to open a stream for writing opr reading
int opencstream(int *perm,char *filename, long int flen){
  int unit=-1;
  int i=0;
  int err=0;
/*   fprintf(stderr, "Length of string %li\n",flen); */
/*   fprintf(stderr, "inputted File name %s\n",filename); */

  //Search for an open spot
  while(unit <1){
    if(i > maxf){
      fprintf(stderr,"opencstream: ERROR too many files open, exiting\n");
      exit(1);
    }
    if(tracker[i] == 0){

      unit=i+1;
      tracker[i]=1;
      
    }else{ 
      i++;
    }
  }

  switch(*perm){
  case 0: // open for  read only
    filelist[unit-1]=fopen(filename,"rb");
    break;
  case 1: //open for write only 
    filelist[unit-1]=fopen(filename,"wb");
    break;
  case 3: //read from standard input
    filelist[unit-1]=stdin;
    break;
  case 4: //write to standard output
    filelist[unit-1]=stdout;
    break;
  }
  if(filelist[unit-1]==NULL){
    fprintf(stderr,"opencstream: ERROR opening file, exiting\n");
    exit(1);
  }
  
/*   /\*Set buffer size to 8192 *\/ */

/*   err=setvbuf(filelist[unit-1],NULL,_IOFBF,8192); */
/*   if(err != 0){ */
/*     fprintf(stderr,"opencstream: ERROR setting buffer, exiting\n"); */
/*     exit(1); */
/*   } */
  return unit;
  
}



//function to close a file
void cclose(int *unit){
  int err;
  //  fprintf(stderr,"unit %i \n",*unit);
  err=fclose(filelist[*unit-1]);
    if(err != 0){
      fprintf(stderr,"ERROR closing file, exiting \n");
      exit(1);
    }

  tracker[*unit-1]=0;
  
}


//Function to write arbitrary data to a file (integer, doubles, real, characters, etc)
// The actual data is the last argument so that it also allowed to pass characters arrays from fortran
void cwrite( int *unit,int *bytes, const void *input){
    if(fwrite(input,1,*bytes,filelist[*unit-1]) != *bytes){
          fprintf(stderr,"ERROR in cwrite: incomplete write\n");
	  exit(1);
    }
  
}

//function which reads in arbitrary data from a file

void cread( int *unit,int *bytes, void *output){
  if(fread(output,1,*bytes,filelist[*unit-1]) != *bytes){
    fprintf(stderr,"ERROR in cread: incomplete read\n");
    exit(1);
  }

}

//Function to skip forward in a file. Will return an error if the file is not seekable
void cskip(int *unit,int *bytes){
  if(fseek(filelist[*unit-1],*bytes,SEEK_CUR) != 0){
    fprintf(stderr,"ERROR in cskip: file not seekable? (pipe/fifo?)\n");
    exit(1);
  }
}

//debug function which prints the memory adress of a pointer to standard error
void printaddress(void *input){
  fprintf(stderr,"%p\n",input);
}


/* /\* Function to copy a pointer adress *\/ */
/* void copyaddress32_(int* out, int* in){ */
/* /\*   fprintf(stdout,"in,out, %d %d \n",*in,*out); *\/ */
/*   *out=*in; */
/* /\*   fprintf(stdout,"in,out, %d %d \n",*in,*out); *\/ */
/* } */

/* /\* Function to copy a pointer adress *\/ */
/* void copyaddress64_(size_t* in,size_t* out, long int flen){ */
/*     fprintf(stdout,"in,out, %p %p \n",in,out); */
/*     out=in; */
/*   fprintf(stdout,"in,out, %p %p \n",in,out); */
/* } */

void cswap(int * bytes, int * sz, unsigned char* data){
  int half=*bytes/2;
  unsigned char dum;
  int shft;
  int i,j;
  for ( i=0; i< *sz ; i++){
    shft=*bytes *i;
     for( j=0; j< half;j++){
       dum=data[shft+j];
       data[shft+j]=data[shft+*bytes-j-1];
       data[shft+*bytes-j-1]=dum;
            } 
  }

}
