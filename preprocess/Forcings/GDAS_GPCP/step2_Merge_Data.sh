#!/bin/bash
#-------------------------------------------------------------


#Year1="1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017"
mkdir -p GDAS_GPCP
Year1="2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016"
for Year in ${Year1};do
Month1="01 02 03 04 05 06 07 08 09 10 11 12"
for Month in ${Month1};do
DAY=1
EDDAY=`cal ${Month} ${Year} | awk 'NF {DAYS = $NF}; END {print DAYS}'`
echo "Last Day in ${Year} ${Month} is $EDDAY"
while [ ${DAY} -le ${EDDAY} ] ; do
if [ ${DAY} -lt 10 ] ; then
DAY2=0${DAY}
else
DAY2=${DAY}
fi
cdo settaxis,${Year}-${Month}-${DAY2},00:00:00,3hour GLDAS_NOAH025_3H.A${Year}${Month}${DAY2}.0000.021.nc4 GLDAS_NOAH025_3H.A${Year}${Month}${DAY2}_00.nc
cdo settaxis,${Year}-${Month}-${DAY2},03:00:00,3hour GLDAS_NOAH025_3H.A${Year}${Month}${DAY2}.0300.021.nc4 GLDAS_NOAH025_3H.A${Year}${Month}${DAY2}_03.nc
cdo settaxis,${Year}-${Month}-${DAY2},06:00:00,3hour GLDAS_NOAH025_3H.A${Year}${Month}${DAY2}.0600.021.nc4 GLDAS_NOAH025_3H.A${Year}${Month}${DAY2}_06.nc
cdo settaxis,${Year}-${Month}-${DAY2},09:00:00,3hour GLDAS_NOAH025_3H.A${Year}${Month}${DAY2}.0900.021.nc4 GLDAS_NOAH025_3H.A${Year}${Month}${DAY2}_09.nc
cdo settaxis,${Year}-${Month}-${DAY2},12:00:00,3hour GLDAS_NOAH025_3H.A${Year}${Month}${DAY2}.1200.021.nc4 GLDAS_NOAH025_3H.A${Year}${Month}${DAY2}_12.nc
cdo settaxis,${Year}-${Month}-${DAY2},15:00:00,3hour GLDAS_NOAH025_3H.A${Year}${Month}${DAY2}.1500.021.nc4 GLDAS_NOAH025_3H.A${Year}${Month}${DAY2}_15.nc
cdo settaxis,${Year}-${Month}-${DAY2},18:00:00,3hour GLDAS_NOAH025_3H.A${Year}${Month}${DAY2}.1800.021.nc4 GLDAS_NOAH025_3H.A${Year}${Month}${DAY2}_18.nc
cdo settaxis,${Year}-${Month}-${DAY2},21:00:00,3hour GLDAS_NOAH025_3H.A${Year}${Month}${DAY2}.2100.021.nc4 GLDAS_NOAH025_3H.A${Year}${Month}${DAY2}_21.nc
DAY=`expr ${DAY} + 1`
done
done #month
cdo mergetime GLDAS_NOAH025_3H.A${Year}${Month}*_*.nc GLDAS_NOAH025_3H.A${Year}${Month}.nc
cdo selname,LWdown_f_tavg GLDAS_NOAH025_3H.A${Year}${Month}.nc GDAS_GPCP/GLDAS_GDAS_3H_LWdown.${Year}${Month}.nc
cdo selname,SWdown_f_tavg GLDAS_NOAH025_3H.A${Year}${Month}.nc GDAS_GPCP/GLDAS_GDAS_3H_SWdown.${Year}${Month}.nc
cdo selname,Tair_f_inst GLDAS_NOAH025_3H.A${Year}${Month}.nc GDAS_GPCP/GLDAS_GDAS_3H_Tair.${Year}${Month}.nc
cdo selname,Qair_f_inst GLDAS_NOAH025_3H.A${Year}${Month}.nc GDAS_GPCP/GLDAS_GDAS_3H_Qair.${Year}${Month}.nc
cdo selname,Psurf_f_inst GLDAS_NOAH025_3H.A${Year}${Month}.nc GDAS_GPCP/GLDAS_GDAS_3H_Psurf.${Year}${Month}.nc
cdo selname,Wind_f_inst GLDAS_NOAH025_3H.A${Year}${Month}.nc GDAS_GPCP/GLDAS_GDAS_3H_Wind.${Year}${Month}.nc
cdo selname,Rainf_f_tavg GLDAS_NOAH025_3H.A${Year}${Month}.nc GDAS_GPCP/GLDAS_GDAS_3H_tot_prcip.${Year}${Month}.nc
rm GLDAS_NOAH025_3H.A${Year}*.nc
done #year
