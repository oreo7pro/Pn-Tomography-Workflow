#include <stdio.h>
#include <math.h>

//--------------------
#include <string.h>
/*next to lines are for mmap()*/
#include <sys/types.h>
#include <sys/mman.h>
/*next two for lseek()*/
#include <sys/types.h>
#include <unistd.h>
//-------------------------------------
double sqrt(),fabs(),atan(),atan2();
 
#include "data14.h"

double sqrt();
//---------------------
#define num 999
#define num2 327680
#define grid_max 2
#define NUMAX 400
//=======add=======
int array_num;
int numax=NUMAX,aa,jj2,i2,jjnew;
int j2,j3,ifany,grid_time,grid_iter,ncell_sdr;
double xx1,yy1,xx2,yy2;
float mesh;
int nnx,nny;
//=================
float xmin=0.0,xmax=360.0,ymin=-90.0,ymax=90.0;
float xmesh=1.0,ymesh=1.0;
float xmesh0=1.0, ymesh0=1.0;
int nx, ny;
int npha;
struct rec *phase;
//#include "data14.h"

double cos(), sin(), exp2(), pow(), amin(), amax(), exp();
int intt();
void *malloc();
/*caddr_t mmap();*/

//-------------------------------------------------------------------------------------
main(argc, argv)
int argc; char **argv;
{
//================================================================================
	int             k,l, ii, jj,len;
	int             iin = 0;
	int		bb,bbb,bbbb;
	double          wei;
	double          res1, res2, res3, res4;
	double			res2slo,res2cos,res2sin;
	double          time;
	int             fdrec, fdsol;
	int             fdevn, fdstn;
	FILE           *fopen();
	double **u, **x, **v;
	int nstat2,nevn2,nsol2;
	double ran1();
	double temp, temp2;
//=============add============
	int kk=0,knum,aj;
	double ksin[num2],kcos[num2],kx[num2],ky[num2];
	int **raysnum;
	int ncell0,nncell;
	int number[num2],number0[num2];
	int xcell0[num],ycell0[num],xcell_temp[num],ycell_temp[num],xcell_sdr[num],ycell_sdr[num];
	double lencell0[num],lencell_temp[num],lencell_sdr[num];
	double xpoint[num],ypoint[num],xpoint0[num],ypoint0[num],xpoint_temp[num],ypoint_temp[num];
	double xpoint_sdr[num],ypoint_sdr[num];
	//============================
	int ncell;
	int xcell[num],ycell[num];
	double lencell[num];
	double azcell[num];
	double length;
	int ntimes,nparams;
	double anorm;
	struct rec *phase2;
//=========================================================================
	struct solution soln;
	double vel,velbias,velerr;
	double ani,anibias,anierr;
	double dir,dirbias,direrr;
	int i, j;
	int count=0;
	int kntzer=0;
	int rmsnm=0;
	double rms=0.0;
	double rav=0.0;
	double ave=0.0; 
	int fdrec2;
	double xloc,yloc;
	float vn=8.1;
	int nxskip=1,nyskip=1;
	int inew;
	float xn,ynn;
	float aaa;
	getarg(argc,argv,"-region","%f/%f/%f/%f",&xmin,&xmax,&ymin,&ymax);
	
	getarg(argc,argv,"-mesh","%f/%f",&xmesh0,&ymesh0);
	getarg(argc,argv,"-numax","%d",&numax);
	fdrec = open("../data/phases.b", 2);
	phase = (struct rec *) mmap(0,lseek(fdrec,0L,SEEK_END),PROT_READ|PROT_WRITE, MAP_SHARED,fdrec,0);
	npha = lseek(fdrec,0L,SEEK_END)/sizeof(struct rec);

//========================================================================================================
	grid_time=1;
	ifany=1;
	for(grid_time=1;ifany>0&&grid_time<grid_max;)
	{
	aa=0;
	xmesh=pow(2,grid_time)*xmesh0/2.0;
	ymesh=pow(2,grid_time)*ymesh0/2.0;
	nx=(xmax-xmin)*xmesh;
	ny=(ymax-ymin)*ymesh;
		 for(jj2=0;jj2<nx*ny;jj2++)
		 {
		 number[jj2]=0;
		 }
		for(i=0,ii=0;i<npha;i++)
		{
	 	if(phase[i].quality<=0.0) continue;
		traceit(phase[i].slat-ymin, phase[i].slon-xmin, phase[i].rlat-ymin, phase[i].rlon-xmin, 
			phase[i].scale,xpoint,ypoint,xcell,ycell,lencell,&ncell,&length);		
		for(j=0;j<ncell;j++) 
		 {
		 jj=xcell[j]+ycell[j]*nx;		
		 number[jj]++;		
		 }
		}
		for(jj=0;jj<nx*ny&&aa==0;jj++)
		 {
	 	 ifany=0;
	 	 if (number[jj]>=numax)
	 	 {
		 fprintf(stderr,"jj number[jj] %d %d \n",jj,number[jj]);//see the first number >numax
		 aa=1;
		 ifany=1;		
		 grid_time++;		
		 continue;
		 }
		}	
	}
fprintf(stderr,"\n grid_time=%d\n",grid_time);
array_num=(xmax-xmin)*pow(2,grid_time*2)*(ymax-ymin)*xmesh0*ymesh0/4;
raysnum=(int**)malloc(sizeof(int*)*(grid_time));
for(i=0;i<grid_time;i++)
raysnum[i]=(int*)malloc(sizeof(int)*array_num);
	//========================================calculate rays number======================
for(i2=0;i2<grid_time;i2++)
{
	
	xmesh=pow(2,i2+1)*xmesh0/2.0;
	ymesh=pow(2,i2+1)*ymesh0/2.0;
	nx=(xmax-xmin)*xmesh;
	ny=(ymax-ymin)*ymesh;
	for(jj=0;jj<nx*ny;jj++)
	{
	number0[jj]=0;
	}
	for(i=0;i<npha;i++)
	{
	 if(phase[i].quality<=0.0) continue;		
	traceit(phase[i].slat-ymin, phase[i].slon-xmin, phase[i].rlat-ymin, phase[i].rlon-xmin, 
			phase[i].scale,xpoint,ypoint,xcell,ycell,lencell,&ncell,&length);		
		for(j=0;j<ncell;j++) 
		 {
		 jj=xcell[j]+ycell[j]*nx;
		 number0[jj]++;			
		 raysnum[i2][jj]=number0[jj];		
		 }
	}
	
}

	mesh=pow(2,grid_time)*xmesh0/2.0;
	nx=(xmax-xmin)*mesh;
	ny=(ymax-ymin)*mesh;

	kk=(xmax-xmin)*(ymax-ymin)*mesh*mesh;
	fprintf(stderr,"grid,xmin,xmax,ymin,ymax,mesh,vn,nx,ny,kk=\n %d %f %f %f %f  %f %f %d %d %d\n",
		grid_time,xmin,xmax,ymin,ymax,mesh,vn,nx,ny,kk);

for(j=0;j<kk;j++)
{
ksin[j]=0;
kcos[j]=0;
}
	fdrec2=open("../data/sol.b",2);
	FILE * fp0508;
	fp0508=fopen("../run/sol.dat","w2");
FILE * fp0718;
	fp0718=fopen("../run/aniadd.txt","w2");

FILE * fp07183;
	fp07183=fopen("../run/aniori.txt","w2");
FILE * fp07184;
	fp07184=fopen("../run/veuni.txt","w2");
FILE * fp07185;
	fp07185=fopen("../run/aniall.txt","w2");
FILE * fp071852;
	fp071852=fopen("../run/veaniall1.dat","w2");
	FILE * fp1;
	fp1=fopen("../output/sin.txt","w2");

	/*loop thru sol.b*/
	for(i=0;read(fdrec2,&soln,sizeof(struct solution));i++)
	{
	if(soln.nu==0) continue;
	//skip every nskip cells
	if((i%nx)%nxskip!=0) continue;
	if((i/nx)%nyskip!=0) continue;
	bb=0;bbb=0;bbbb=0;
	count++;
	vel=soln.ds;
	fprintf(fp071852,"%d %9.4f %7.7f %7.7f\n",i,soln.ds,soln.sin,soln.cos);
			
	//fprintf(stderr,"vel vn=%f %f\n",vel,vn);
	if(vel==0.0) kntzer++;
		else
		{
		ave += vel;
		rms += vel*vel;
		rav += fabs(vel);
		rmsnm++;
		}
//---------------------------------------------------------------------------------------------------------------------------
			 	vel *= -vn*vn;
				velbias = -vn*vn*(soln.sum-soln.ds);
				velerr = vn*vn*sqrt(soln.sum2);

				//fprintf(stdout,"sum2,sum2cos,sum2sin= %g %g %g\n",soln.sum2,soln.sum2cos,soln.sum2sin);
		
				//ani is the anisotropy velocity perturbation
				ani= vn*vn*sqrt(soln.cos*soln.cos+soln.sin*soln.sin);
				anibias = vn*vn*sqrt(soln.sumcos*soln.sumcos + soln.sumsin*soln.sumsin);
				anibias -= ani;
				anierr = soln.sum2cos*soln.cos*soln.cos + soln.sum2sin*soln.sin*soln.sin + 2*soln.sum2cossin*soln.cos*soln.sin;
				anierr = vn*vn*sqrt(anierr)/sqrt(soln.cos*soln.cos+soln.sin*soln.sin);
		if(sqrt(soln.cos*soln.cos+soln.sin*soln.sin)==0)
		fprintf(fp1,"%d %f %f %f\n",i,soln.cos,soln.sin,sqrt(soln.cos*soln.cos+soln.sin*soln.sin));		
				//dir is the atan(sin/cos) w/ coefficient 180/pi/2
				dir=(28.6479)*atan2(soln.sin,(soln.cos+.000000000001));
				dir = 90.0 - dir; //convert to gmt coords (pos counterclockwise from E
				dir -= 90.0; //plot max velocity not max slowness
//======================================




				dirbias = (28.6479)*atan2(soln.sumsin,soln.sumcos+.000000000001);
				dirbias = 90.0 - dirbias; //convert to gmt coords (pos counterclockwise from E
				dirbias -= 90.0; //plot max velocity not max slowness
				dirbias -= dir;
				direrr = 28.6419*sqrt(soln.sum2sin*soln.cos*soln.cos + soln.sum2cos*soln.sin*soln.sin + 2*soln.sum2cossin*soln.cos*soln.sin);
				direrr /= (soln.cos*soln.cos + soln.sin*soln.sin);
//============================================zero ksin kcos==================================

	//fprintf(fp0508,"%d %f\n",i,vel);
	j=0;	
	aa=pow(2,grid_time-1);
	inew=(int)(   (float)(i%nx)/aa+(int)((float)(i/nx)/aa)  *nx/aa    );
	if(  (i%nx)%aa==0 && (i/nx)%aa==0 && raysnum[j][inew]<numax  )
		{
		
		//aa=aa*2;
//---------------------------------------------------------------------------------------------------
			  for(xn=0;xn<aa;xn++)
			  {
			  xloc=xmin + (float)(i%nx)/mesh +(float)0.5/mesh+ (float)(1*xn)/mesh;
				for(ynn=0;ynn<aa;ynn++)
			  	{
				yloc=ymin + (float)(i/nx)/mesh+(float)0.5/mesh + (float)(1*ynn)/mesh;
				fprintf(fp0508,"%9.4f %9.4f %5d %7.7f %7.7f  %7.7f %7.7f %7.7f %7.7f %7.7f  %7.7f %7.7f\n",
				xloc,yloc,soln.nu,vel,velbias,
				velerr,ani,anibias,anierr,dir,
				dirbias,direrr);
				}
			  }
//==================================================================anisotropy
			xloc=xmin + (float)(i%nx)/mesh +(float)1/mesh;
			yloc=ymin + (float)(i/nx)/mesh +(float)1/mesh;
			fprintf(fp07183,"%9.4f %9.4f %5d %7.7f %7.7f\n",xloc,yloc,soln.nu,ani,dir);
			fprintf(fp07185,"%9.4f %9.4f %5d %7.7f %7.7f\n",xloc,yloc,soln.nu,ani,dir);
			
			fprintf(fp07184,"%9.4f %9.4f %5d %7.7f \n",xloc,yloc,soln.nu,vel);
		}
		
	if(raysnum[j][inew]>=numax)
		{
		xloc=xmin + (float)(i%nx)/mesh + (float)0.5/mesh;
		yloc=ymin + (float)(i/nx)/mesh + (float)0.5/mesh;
		
		
		fprintf(fp0508,"%9.4f %9.4f %5d %7.7f %7.7f  %7.7f %7.7f %7.7f %7.7f %7.7f  %7.7f %7.7f\n",
				xloc,yloc,soln.nu,vel,velbias,
				velerr,ani,anibias,anierr,dir,
				dirbias,direrr);
		fprintf(fp07183,"%9.4f %9.4f %5d %7.7f %7.7f\n",xloc,yloc,soln.nu,ani,dir);
		if(i%2==0&&i/nx%2==0)
		fprintf(fp07185,"%9.4f %9.4f %5d %7.7f %7.7f\n",xloc+ (float)0.5/mesh,yloc+ (float)0.5/mesh,soln.nu,ani,dir);
		
fprintf(fp07184,"%9.4f %9.4f %5d %7.7f \n",xloc,yloc,soln.nu,vel);
		knum=(i%nx)/2+((i/nx)/2)*nx/2;
		ksin[knum] = ksin[knum]+ ani*sin(dir*3.1415926/180);
		kcos[knum] = kcos[knum]+ ani*cos(dir*3.1415926/180);
		}

}
fdrec2=open("../data/sol.b",2);
/*loop thru sol.b*/
for(i=0;read(fdrec2,&soln,sizeof(struct solution));i++)
{
if(soln.nu==0) continue;
	
if((i%nx)%nxskip!=0) continue;
if((i/nx)%nyskip!=0) continue;

count++;		
//----------------------------------------------
j=0;
aa=pow(2,grid_time-1);
inew=(int)(   (float)(i%nx)/aa+(int)((float)(i/nx)/aa)  *nx/aa    );
if(  raysnum[j][inew]>=numax  )
{
	
		//knum=(i%nx)/2+(i/nx)*(nx/4);
		knum=(i%nx)/2+((i/nx)/2)*nx/2;
		xloc=xmin + (float)(knum%(nx/2))*2/mesh + (float)0.5*2/mesh;
		yloc=ymin + (float)(knum/(nx/2))*2/mesh + (float)0.5*2/mesh;
		ani=(float)(sqrt(ksin[knum]*ksin[knum]+kcos[knum]*kcos[knum]))/4;
		dir=57.2958*atan2(ksin[knum],kcos[knum]);		
		//fprintf(fp07182,"%9.4f %9.4f %5d %7.7f %7.7f\n",xloc,yloc,soln.nu,ani,dir);

	}
}

	for(ii=0;ii<grid_time;ii++)
	free(raysnum[ii]);
	free(raysnum);
	fclose(fp1);
	fclose(fp0508);
	fclose(fp0718);
	fclose(fp07183);
	fclose(fp07184);
	fclose(fp07185);
	fclose(fp071852);
	fprintf(stderr,"\n %d zeros of %d total\n",kntzer,count);
	fprintf(stderr,"i= %d\n",i);
	fprintf(stderr,"rms= %f , rmsnm= %d \n",sqrt(rms/((float)rmsnm+.0001)),rmsnm);
	fprintf(stderr,"aveabs= %f , rmsnm= %d \n",(rav/((float)rmsnm+.0001)),rmsnm);
	fprintf(stderr,"ave = %f\n\n",ave/(rmsnm+.0001));
}
//------------------------------------------------------------------------------------------------


traceit(slat, slon, rlat, rlon, scale,xpoint,ypoint, xcell,ycell,lencell,ncell,length)
/* computes the travel time from the locations */
/* assumes straight lines */
/*xcell,ycell,are matrixes that contain the cell number*/
/**ncell is the total number of cells the ray crosses*/
float	slat, slon, rlat, rlon, scale;
int *ncell;
int *xcell,*ycell;
double *xpoint,*ypoint;
double *lencell, *length;
{	
	double	x1, y1, x2, y2, dx, dy;
	double	len, x, y;
	int		ix, iy;
	double	xnew, ynew;
	int		ixnew, iynew;
	double	alpha, sgnalp;
	double	amax();
	
	/*zero some stuff*/
	*ncell=0;
	
	/* set rlon to be minimum, if not switch */
	if(slon<rlon) {
		x1 = x = slon;
		x2 = rlon;
		y1 = y = slat;
		y2 = rlat;
	} else {
		x1 = x = rlon;
		x2 = slon;
		y1 = y = rlat;
		y2 = slat;
	}

	/* set some parameters */
	dx = x2 - x1;
	dy = y2 - y1;
	if(dy == 0.0) dy = 0.0001;
	*length = sqrt(dx*dx + dy*dy) * scale;
	ix = intt(x1 * xmesh);
	iy = intt(y1 * ymesh);
	alpha = dx / dy;
	sgnalp = sign(alpha);
	
	//fprintf(stderr,"x,x2=%f %f\n",x,x2);
	/* now trace */
	while (x<x2) {
		ynew = (iy + amax(sgnalp, 0.0)) / ymesh;
		xnew = x + alpha * (ynew - y);
		if(xnew>x2) {
			xnew = x2;
			ynew = y2;
		}
		ixnew = intt(xnew * xmesh);
		iynew = iy + sgnalp;
		if(ixnew>ix) {
			ixnew = ix + 1;
			xnew = ((float) ixnew) / xmesh;
			ynew = y + (xnew - x) / alpha;
			iynew = iy;
		}
		dx = xnew - x;
		dy = ynew - y;
		len = sqrt(dx*dx + dy*dy) * scale;
		/* this is the important loop */
		if(ix >= 0 && ix<nx && iy >= 0 && iy<ny) {
			xcell[*ncell]=ix;
			ycell[*ncell]=iy;
			xpoint[*ncell]=x;
			ypoint[*ncell]=y;
			lencell[*ncell]=len;
			(*ncell)++;
		}
		
		ix = ixnew;
		iy = iynew;
		x = xnew;
		y = ynew;
		xpoint[*ncell]=x;
		ypoint[*ncell]=y;
	}
	return;
}


