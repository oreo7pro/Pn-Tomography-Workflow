#此文件用于绘制带sks的各向异性图件

#!sh
set -x
#has one option which is a file containing the parameters amd a file name

#first read in the important parameters from the optfile

REGION=15/53/25/48
RVALUE=20/48/30/43
JVALUE=m0.8

CELLSIZE=15m/15m
AVEVEL=8.05

ps=sks各向异性.ps
FILEDAT2=aniall.txt
#FILEDAT2=aniori.txt
FILEDAT3=output_avg0.5.dat
COLORFILE=color30.cpt
FILE=../data/sol.b


gmt grdimage sol_clip.grd -C$COLORFILE -J$JVALUE -R$RVALUE -V -K  -X3.5 -Y4.5  -Ba5f2.5 >$ps
#横版全图
#gmt grdimage sol_clip.grd -C$COLORFILE -J$JVALUE -R$RVALUE -V -K  -X2.8 -Y5.8  -Ba1f1 -P >$ps
#-P竖过来 竖版小图
#-Itopo_grad_eq2.grd -Cstation.cpt
 

gmt pscoast -A1000/1 -Di  -N1/0.35p -N2/0.25p -J$JVALUE  -Rg$RVALUE -W0.75p,white   -K -O>>$ps

#awk '{if($3>0) print $1" "$2" "$5" "$4/1.8}' $FILEDAT2 |gmt psxy -L -W0.75p  -J$JVALUE -R$RVALUE -Sv1c -V -K -O>>$ps
#awk '{if($3>0) print $1" "$2" "$5+180" "$4/1.8}' $FILEDAT2 |gmt psxy -L -W0.75p -J$JVALUE -R$RVALUE -Sv1c -V -K -O >>$ps
awk '{if($3>0) print $1" "$2" "$5" "$4/2.3}' $FILEDAT2 |gmt psxy -L -W0.25p  -J$JVALUE -R$RVALUE -Sv1c -V -K -O>>$ps
awk '{if($3>0) print $1" "$2" "$5+180" "$4/2.3}' $FILEDAT2 |gmt psxy -L -W0.25p -J$JVALUE -R$RVALUE -Sv1c -V -K -O >>$ps
##########################################
awk '{if($4>=0.8&&$4<1.1) print $1" "$2" "$3" "$4/9.6}' $FILEDAT3 |gmt psxy  -L -J$JVALUE -R$RVALUE -SV1c -W2p -V -K -O >>$ps
awk '{if($4>=0.8&&$4<1.1) print $1" "$2" "$3+180" "$4/9.6}' $FILEDAT3 |gmt psxy  -L -J$JVALUE -R$RVALUE -SV1c -W2p -V -K -O >>$ps

awk '{if($4>=0.8&&$4<1.1) print $1" "$2" "$3" "$4/11.34}' $FILEDAT3 |gmt psxy  -L -J$JVALUE -R$RVALUE -SV1c -W1p,white -V -K -O >>$ps
awk '{if($4>=0.8&&$4<1.1) print $1" "$2" "$3+180" "$4/11.34}' $FILEDAT3 |gmt psxy  -L -J$JVALUE -R$RVALUE -SV1c -W1p,white -V -K -O >>$ps


############################
awk '{if($4>=1.1&&$4<1.4) print $1" "$2" "$3" "$4/9.6}' $FILEDAT3 |gmt psxy  -L -J$JVALUE -R$RVALUE -SV1c -W2p -V -K -O >>$ps
awk '{if($4>=1.1&&$4<1.4) print $1" "$2" "$3+180" "$4/9.6}' $FILEDAT3 |gmt psxy  -L -J$JVALUE -R$RVALUE -SV1c -W2p -V -K -O >>$ps

awk '{if($4>=1.1&&$4<1.4) print $1" "$2" "$3" "$4/10.99}' $FILEDAT3 |gmt psxy  -L -J$JVALUE -R$RVALUE -SV1c -W1p,white -V -K -O >>$ps
awk '{if($4>=1.1&&$4<1.4) print $1" "$2" "$3+180" "$4/10.99}' $FILEDAT3 |gmt psxy  -L -J$JVALUE -R$RVALUE -SV1c -W1p,white -V -K -O >>$ps
############################

awk '{if($4>=1.4&&$4<1.7) print $1" "$2" "$3" "$4/9.6}' $FILEDAT3 |gmt psxy  -L -J$JVALUE -R$RVALUE -SV1c -W2p -V -K -O >>$ps
awk '{if($4>=1.4&&$4<1.7) print $1" "$2" "$3+180" "$4/9.6}' $FILEDAT3 |gmt psxy  -L -J$JVALUE -R$RVALUE -SV1c -W2p -V -K -O >>$ps

awk '{if($4>=1.4&&$4<1.7) print $1" "$2" "$3" "$4/10.68}' $FILEDAT3 |gmt psxy  -L -J$JVALUE -R$RVALUE -SV1c -W1p,white -V -K -O >>$ps
awk '{if($4>=1.4&&$4<1.7) print $1" "$2" "$3+180" "$4/10.68}' $FILEDAT3 |gmt psxy  -L -J$JVALUE -R$RVALUE -SV1c -W1p,white -V -K -O >>$ps
###########################
awk '{if($4>=1.7) print $1" "$2" "$3" "$4/9.6}' $FILEDAT3 |gmt psxy  -L -J$JVALUE -R$RVALUE -SV1c -W2p -V -K -O >>$ps
awk '{if($4>=1.7) print $1" "$2" "$3+180" "$4/9.6}' $FILEDAT3 |gmt psxy  -L -J$JVALUE -R$RVALUE -SV1c -W2p -V -K -O >>$ps

awk '{if($4>=1.7) print $1" "$2" "$3" "$4/10.26}' $FILEDAT3 |gmt psxy  -L -J$JVALUE -R$RVALUE -SV1c -W1p,white -V -K -O >>$ps
awk '{if($4>=1.7) print $1" "$2" "$3+180" "$4/10.26}' $FILEDAT3 |gmt psxy  -L -J$JVALUE -R$RVALUE -SV1c -W1p,white -V -K -O >>$ps

awk '{if($4>=0.8&&$4<1.5) print $1" "$2}' $FILEDAT3 |gmt psxy -L -J$JVALUE -R$RVALUE -Gpink -W0.004c -Sc0.06c -V -K -O >> $ps
awk '{if($4>=1.5) print $1" "$2}' $FILEDAT3 |gmt psxy -L -J$JVALUE -R$RVALUE -Gpink -W0.004c -Sc0.066c -V -K -O >> $ps



gmt psxy  -L -J$JVALUE -R$RVALUE -SV1c -W2p  <<END -O -K >>$ps
128 20 45 0.104167
128 20 225 0.104167
128.5 20 45 0.15625
128.5 20 225 0.15625
129 20 45 0.208333
129 20 225 0.208333
END
gmt psxy  -L -J$JVALUE -R$RVALUE -SV1c -W1p,white  <<END -V -K -O >>$ps
128 20 45 0.08818
128 20 225 0.08818
128.5 20 45 0.140449
128.5 20 225 0.140449
129 20 45 0.19493
129 20 225 0.19493
END

gmt psxy -L -J$JVALUE -R$RVALUE -Gpink -W0.004c -Sc0.06c <<END  -V -K -O >> $ps
128 20 
128.5 20
129 20
END

#awk '{if($3>0) print $1" "$2" "$5" "$4/2.7}' aniori.txt |gmt psxy -L -W0.05p -J$JVALUE -R$RVALUE -Sv1c -V -K -O>>$ps

#awk '{if($3>0) print $1" "$2" "$5+180" "$4/2.7}' aniori.txt |gmt psxy -L -W0.05p -J$JVALUE -R$RVALUE -Sv1c -V -K -O >>$ps
#first get the data and grid it

gmt pscoast -A1000/1 -Di  -J$JVALUE -R$RVALUE -W0.5p  -K -O>>$ps

awk '{ print $1" "$2}' world_volcanoes.dat |gmt psxy  -J$JVALUE -R$RVALUE -V -Gred -W0.3p -St0.3c -K -O  >> $ps

gmt psxy  -Jx1 -R0/20/0/10 -V -Sv0.4c -W1p -L  <<END -O -K -X-0.25 -Y-5>>$ps
5.0 1.5 45 0.00
7.25 1.5 45 0.09
8.9 1.5 45 0.18
10.3 1.5 45 0.27 
END
gmt pstext -Jx1 -R0/25/0/10 -F+f20p -Y-0.3  -V <<END -O >>$ps
5.35     .96   @!0
5.95    .96   @!%

6.95     .96   @!1
7.48     .96   @!%

8.5     .96   @!2
9.1     .96   @!%

10.05     .96   @!3
10.65     .96   @!%

8.0 0.5 Pn Anisotropy
END


###########################
echo "all done now"

