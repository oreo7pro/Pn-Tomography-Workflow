#!/bin/sh
#this shell is for running ctomo
#	it copies from sta0.b etc
#	it does not do test models
set -x
date

#first read in the important parameters from the optfile

#East_Med
REGION=10/58/20/53

MESH=4/4
date

#now invert15
../bin/ctomo -region $REGION -mesh $MESH -nlsqr 200 -weight 500 -aniwei 500 -numax 800 -nbootv 0
date

exit