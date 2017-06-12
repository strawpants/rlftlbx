#!/bin/bash

fdeg1=Swensondeg1.txt

wget -O $fdeg1 ftp://podaac.jpl.nasa.gov/allData/tellus/L2/degree_1/deg1_coef.txt

#convert degree 1 potential coefficients from
sc=11047256.4063275

sed -r -e 's/[0-9]\s+1\s+[01]/&\ /g' -e '/^[^2]/ d' -e '/^$/ d' -e 's/-0\./\ -0\./g' $fdeg1  | awk '{yr=substr($1,1,4);mn=substr($1,5,2)}/[0-9] +1 +0/{z='$sc'*$4;sz='$sc'*$6}/[0-9] +1 +1/{x='$sc'*$4;y='$sc'*$5;sx='$sc'*$6;sy='$sc'*$7;printf "15-%02d-%4d %e %e %e %e %e %e\n",mn,yr,x,y,z,sx,sy,sz}' | GPS_calendar -y -D=1/10 | sed -e 's/^[ \t]*//' > SwensonCM-CF.txt
