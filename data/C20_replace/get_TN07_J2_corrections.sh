#!/bin/bash

#bash script to convert the Ries technical note to something readable by c20_replace.f90

#note: the time tags are shifted by 15 days

#updated 27 Jan 2015 (include download link)

fout=C20_TNries.txt
#backup old file
cp $fout ${fout}.bk

#get new file online
ftplink=ftp://podaac.jpl.nasa.gov/allData/grace/docs/TN-07_C20_SLR.txt
wget -N $ftplink


fin=`basename $ftplink`

awk 'BEGIN{p=0}/^PRODUCT/{p=1;l=NR}p==1&&NR>l+1{printf "echo `GPS_calendar -y \x27 MJD %f\x27` %e %e\n", $1+15,$4*1e-10,$5*1e-10}' $fin | bash > $fout
