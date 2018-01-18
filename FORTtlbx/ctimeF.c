/*Fortran interface for ctime.h*/
/*Author Roelof Rietbroek*/
/*date 18 Jan 2018*/

#include <time.h>
#include <stdio.h>
/*Returns a C style tm struct which contains date and time */
time_t fmktime_(int * yy, int * mm, int * dd, int * hh, int * min, int * ss){
    struct tm datetime;
    tzset();//initialize current time zone
    /*fprintf(stderr,"%d\n",timezone);*/
    //note we correct for the local time zone in order to work with GMT times
    datetime.tm_sec=*ss-timezone;  //int seconds after the minute    0-61
    datetime.tm_min=*min;  //int minutes after the hour  0-59
    datetime.tm_hour=*hh; //int hours since midnight    0-23
    datetime.tm_mday=*dd; //int day of the month    1-31
    datetime.tm_mon=*mm; //  int months since January    0-11
    datetime.tm_year=*yy-1900; //int years since 1900    
    
    //the following 2 values are not needed (will be recomputed in mktime)
    //datetime.tm_wday=0; //int days since Sunday   0-6
    //datetime.tm_yday=0; //int days since January 1    0-365
    
    
    datetime.tm_isdst=0; //int Daylight Saving Time flag   
    
    //note mktime assume the local time so that's wgy we corrected above with  
    time_t datet=mktime(&datetime);

    //assume GMT (correct for current timezone)
    return datet;
} 

//function to format a time in a string
void fstrftime_(char* outputstr,char* frmt, time_t * datetime,int flen1, int flen2){
    /*fprintf(stderr,"input lengths %d %d\n",flen1,flen2);*/
    strftime (outputstr, flen2, frmt, localtime(datetime) );    


}




