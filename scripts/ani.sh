#!sh
set -x
#this plots the velocity map
#has one option which is a file containing the parameters amd a file name

#first read in the important parameters from the optfile

#East_Med
REGION=10/58/20/53
RVALUE=20/48/30/43
JVALUE=m0.8

# === 新增：棋盘格网格线控制参数 ===
# 通过调整 LON_START 和 LAT_START，你可以随意平移整个网格
# 如果你的棋盘格边界正好是从整数经纬度开始的，保持 20 和 30 即可
GRID_SPACING=1.25
LON_START=20.12
LAT_START=30.12

CELLSIZE=15m/15m
AVEVEL=8.05
ps=Ani_noise10.ps
FILEDAT3=aniall_noise10.txt

COLORFILE=color30.cpt
FILE=../data/sol.b

gmt grdimage sol_clip.grd -C$COLORFILE -J$JVALUE -R$RVALUE -V -K  -X4.5 -Y4.8  -Ba5f2.5 >$ps
# === 新增：生成并绘制 1.25° 黑色细线网格 ===
# 提取 RVALUE 中的经纬度边界
MIN_LON=20
MAX_LON=48
MIN_LAT=30
MAX_LAT=43

# 使用 awk 自动生成经纬度网格线坐标
awk -v lon0=$LON_START -v lat0=$LAT_START -v step=$GRID_SPACING \
    -v min_lon=$MIN_LON -v max_lon=$MAX_LON -v min_lat=$MIN_LAT -v max_lat=$MAX_LAT 'BEGIN {
    # 生成经度线 (竖线)
    for(lon=lon0; lon<=max_lon; lon+=step) {
        print lon, min_lat
        print lon, max_lat
        print ">"
    }
    # 生成纬度线 (横线)
    for(lat=lat0; lat<=max_lat; lat+=step) {
        print min_lon, lat
        print max_lon, lat
        print ">"
    }
}' > checker_grid.txt

# 用 psxy 将网格线画在底图上 (W0.25p,black 表示0.25磅黑色细线)
#gmt psxy checker_grid.txt -J$JVALUE -R$RVALUE -W0.25p,black -A -V -K -O >> $ps
rm -f checker_grid.txt
# ===========================================


#-Itopo_grad_eq2.grd -Cstation.cpt
awk '{if($3>0) print $1" "$2" "$5" "$4/2.5}' $FILEDAT3 |gmt psxy -L -W0.25p  -J$JVALUE -R$RVALUE -Sv1c -V -K -O>>$ps
awk '{if($3>0) print $1" "$2" "$5+180" "$4/2.5}' $FILEDAT3 |gmt psxy -L -W0.25p -J$JVALUE -R$RVALUE -Sv1c -V -K -O >>$ps


gmt pscoast -A1000/1 -Di -N1/0.3p -J$JVALUE -R$RVALUE -W0.45p  -K -O>>$ps
#awk '{ print $1" "$2}' world_volcanoes.dat |gmt psxy  -J$JVALUE -R$RVALUE -V -Gred -W0.3p -St0.3c -K -O  >> $ps

gmt psxy  -Jx1 -R0/20/0/10 -V -Sv0.4c -W1p -L  <<END -O -K -X1 -Y-3.5>>$ps
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

