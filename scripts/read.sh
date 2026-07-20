set -x

#East_Med
REGION=10/58/20/53


MESH=4/4
FILE=../data/sol.b
NUMAX=800

../bin/readsol -file $FILE -region $REGION -mesh $MESH -numax $NUMAX > nonuni2.dat
#../bin/readsol -file $FILE -region $REGION -mesh $MESH > sol.dat
echo "all done now"
