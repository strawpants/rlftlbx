/*C methods which handle regular expressions queries
These routines are constructed such that they can be called from Fortran through a fortran wrapper around C.
The wrapper is important as it puts null characters at the appropriate places in the strings
*/

/*Coded by Roelof Rietbroek 02-03-2009*/


#include <regex.h>
#include <stdio.h>
#include <stdlib.h>

/*Declare a list of 300 compiled regular exppressions */
#define MAXR 2000
/*const int maxr=2000;*/
regex_t reglist[MAXR]; //Maximum of maxr compiled regular expressions
int regtracker[MAXR] = {0}; //tracker which keeps hold of the compiled regular expressions

/* Function which compiles a new regular expression and puts it in an empty slot */

int compregex_(char *regex, long int flen){
  int unit=-1;
  int i=0;
  int err=0;
  //Search for an open spot
  while(unit <1){
    if(i > MAXR){
      fprintf(stderr,"compregex: ERROR too many compiled regular expressions, exiting\n");
      exit(1);
    }
    if(regtracker[i] == 0){

      unit=i+1;
      regtracker[i]=1;
      
    }else{ 
      i++;
    }
  }


  /*Compile regular expression  */
  
  err=regcomp(&reglist[unit-1], regex, REG_EXTENDED);
  if(err!=0){
    fprintf(stderr,"compregex: ERROR compiling regular expression, exiting\n");
    exit(1);
  }
/*   fprintf(stderr,"compregex: %i %s \n",unit,regex); */
  return unit; // return the entry of the compiled expression
}

int execregex_(int *regnum, char *strng, long int flen){
 /*     fprintf(stderr,"compregex: %s \n",strng); */
  /* return 0 if the string matched the regular expression or 1 other wise */
  return regexec(&reglist[*regnum-1], strng, (size_t) 0, NULL, 0);

}
