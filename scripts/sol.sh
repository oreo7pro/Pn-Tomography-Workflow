#!sh
set -x
#this plots the velocity map
#has one option which is a file containing the parameters amd a file name

#first read in the important parameters from the optfile

#East_Med
REGION=10/58/20/53
RVALUE=20/48/30/43
JVALUE=m0.8

CELLSIZE=15m/15m

FILEDAT=sol_noise10.dat

ps=Velo_noise10.ps

COLORFILE=color30.cpt

#gmt makecpt -Cpolar -T-0.6/0.6 >custom.cpt

awk '{ if($3>0) print $1,$2,$4}' $FILEDAT|gmt xyz2grd -Am -r -Gsol.grd -I$CELLSIZE -R$REGION -V
gmt grdinfo sol.grd -L2

#切换到新区域的时候重新运行一遍这部分
#gmt grdgradient /home/heyuhui/Desktop/nonuiform/run/earth_relief_15s.grd -A0.0 -R$RVALUE -Gtopo_grad.grd -V
#gmt grdhisteq topo_grad.grd -Gtopo_grad_eq1.grd -N
#gmt grdmath topo_grad_eq1.grd 0.25 MUL = topo_grad_eq1.grd


#clip
# -Itopo_grad.grd 
gmt grdsample sol.grd -Gresol.grd -I1m
gmt grdsample topo_grad_eq1.grd -Gtopo_grad_eq2.grd -I1m
#grdclip resol.grd -Gsol_clip.grd -Sa8.55/8.55 -Sb7.55/7.55
gmt grdclip resol.grd -Gsol_clip.grd -Sa0.5/0.5 -Sb-0.5/-0.5

#gmt grdimage sol_clip.grd    -C$COLORFILE -J$JVALUE -R$RVALUE -V    -X3.5 -Y4.5 -Ba5f2.5 -P -K  >$ps
gmt grdimage sol_clip.grd -Itopo_grad_eq2.grd -C$COLORFILE -J$JVALUE -R$RVALUE -V  -X4.5 -Y4.8 -Ba5f2.5  -K  >$ps

awk '{ print $1" "$2}' world_volcanoes.dat |gmt psxy  -J$JVALUE -R$RVALUE -V -Gred -W0.3p -St0.3c -K -O  >> $ps 

#gmt psxy PB2002_boundaries2.gmt -J$JVALUE -R$RVALUE -V -W1p,54/54/54  -K -O  >> $ps

gmt pscoast -A1000/1 -Di  -N1/0.3p -J$JVALUE  -Rg$RVALUE -W0.4p   -K -O>>$ps


#this part to draw South China Sea insert
if [ ${J2} ]
	then
	sh scs.sh sol.ps
fi

gmt psscale  -C$COLORFILE -Bxa0.125 -D10.87/-1.5/9.0/0.5h -K  -O -X-1.6 -Y-0>>$ps
gmt pstext  -Jx1 -R0/15/0/8.4  -Y-4.5 -V <<END -O  >>$ps
10.85 1.2 20 0.0 1 10  Pn Anomaly (km/s) 
END


echo "all done now"

