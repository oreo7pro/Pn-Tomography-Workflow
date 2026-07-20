#!/bin/sh 
#this shell is for running rdphase
#ctomo.sh is jst to use rdphase
set -x
date
#first read in the important parameters from the optfile

#East_Med
REGION=10/58/20/53

#Aegean
#REGION=15/35/27/47

SCR=../rawdata
INFILES="$SCR/Pn6420_Eurasia_g_cor50_m50_new2021.dat"

RDPHASE="-stnmin 5 -evnmin 5 -resmax 5 -cresmax 5 -delmax 12 -delmin 2 -depmax 50"
#-resmax 5 

#next get the data set together 
cat $INFILES | ../bin/rdphase -region $REGION $RDPHASE
exit

