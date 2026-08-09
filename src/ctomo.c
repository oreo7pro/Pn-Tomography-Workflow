/* ctomo.c - to invert travel times , this version also updates the dtimes in phases.b */
/* this version for Pn, it accounts for the ray offset */
/* this does gauss siedel iteration */
/* non-uniform grids */
/* & this version does not do the scaleing */
/* & this version bins the stations */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
/*next to lines are for mmap()*/
#include <sys/types.h>
#include <sys/mman.h>
/*next two for lseek()*/
#include <sys/types.h>
#include <unistd.h>



/* inversion parameters*/
#define NLSQR 200			/*number of lsqr iterations*/
#define NBOOTV	0			/*number of bootstrap iterations for variance*/
#define WEIGHT	400		/*amplitude of laplacian damping*/
#define ANIWEI	400		/*include anisotropy, weight for anisotropy*/
#define PIO180	0.01745329252
#define pi 3.1415926
#define N 999       
#define Rr 6371
#define num 999
#define num2 327680
#define grid_max 2
#define NUMAX 800
//=======add=======
int array_num;
int mesh_temp,aa,jj2,i2,jjnew;
int j2,j3,ifany,grid_time,grid_iter,ncell_sdr;
double xx1,yy1,xx2,yy2;
float mesh;
int nnx,nny;
//=================
float xmin=0.0,xmax=360.0,ymin=-90.0,ymax=90.0;
float xmesh=1.0,ymesh=1.0;
float xmesh0=1.0,ymesh0=1.0;
float weight=WEIGHT, aniwei=ANIWEI;
int numax=NUMAX;
int nlsqr=NLSQR, nbootv=NBOOTV;
int nx, ny;
int npha;
int nstat, nevn;
struct rec *phase;
double *utemp,*xtemp,*vtemp;

int **raysnum;

#include "data14.h"

double sqrt(), fabs(), cos(), sin(), exp2(), pow(), amin(), amax(), exp();
int intt();
void *malloc();
/*caddr_t mmap();*/


struct locations *station, *event;
struct solution *sol;

main(argc,argv)
int argc; char **argv;
{
	int             i, j, k,l, ii, jj,len;
	int             iin = 0;
	
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
	
	int ncell0,nncell;
	int number[num2],number0[num2];
	int xcell0[num],ycell0[num],xcell_temp[num],ycell_temp[num],xcell_sdr[num],ycell_sdr[num];
	double lencell0[num],lencell_temp[num],lencell_sdr[num];
	double xpoint[num],ypoint[num],xpoint0[num],ypoint0[num],xpoint_temp[num],ypoint_temp[num];
	double xpoint_sdr[num],ypoint_sdr[num];
	int jj1,jj11,jj2,jj22;
	//============================
	int ncell;
	int xcell[num],ycell[num];
	double lencell[num];
	double azcell[num];
	double length;
	int ntimes,nparams;
	double anorm;
	struct rec *phase2;
	
	if(getarg(argc,argv,"-region","%f/%f/%f/%f",&xmin,&xmax,&ymin,&ymax)<4) errprint("ctomo: cant getarg region");
	getarg(argc,argv,"-weight","%f",&weight); 
	getarg(argc,argv,"-aniwei","%f",&aniwei); 
	getarg(argc,argv,"-numax","%d",&numax);
	getarg(argc,argv,"-nlsqr","%d",&nlsqr);
	getarg(argc,argv,"-nbootv","%d",&nbootv);
	getarg(argc,argv,"-mesh","%f/%f",&xmesh0,&ymesh0);
	if(nbootv!=NBOOTV&&NBOOTV*nbootv==0) errprint ("ctomo: must recompile with NBOOTV reset");
	if(aniwei!=ANIWEI&&ANIWEI*aniwei==0) errprint ("ctomo: must recompile with ANIWEI reset");
	fprintf(stderr,"xmin,xmax,ymin,ymax,xmesh0,ymesh0= %f %f %f %f %f %f\n",xmin,xmax,ymin,ymax,xmesh0,ymesh0);
	fprintf(stderr,"weight,aniwei,nlsqr,nbootv= %f %f %d %d\n",weight,aniwei,nlsqr,nbootv);

	
	/* open all files*/
	fdrec = open("../data/phases.b", 2);
	fdevn = open("../data/events.b", 2);
	fdstn = open("../data/stations.b", 2);
	fdsol = creat("../data/sol.b", 0644);
	
	
	fprintf(stderr,"fdrec,fdevn,fdstn,fdsol= %d %d %d %d\n", fdrec, fdevn, fdstn, fdsol);

	/* mmap memory spaces */
	phase = (struct rec *) mmap(0,lseek(fdrec,0L,SEEK_END),PROT_READ|PROT_WRITE, MAP_SHARED,fdrec,0);
	npha = lseek(fdrec,0L,SEEK_END)/sizeof(struct rec);
	station = (struct locations *) mmap(0,lseek(fdstn,0L,SEEK_END),PROT_READ|PROT_WRITE, MAP_SHARED,fdstn,0);
	nstat = lseek(fdstn,0L,SEEK_END)/sizeof(struct locations);
	event = (struct locations *) mmap(0,lseek(fdevn,0L,SEEK_END),PROT_READ|PROT_WRITE, MAP_SHARED,fdevn,0);
	nevn = lseek(fdevn,0L,SEEK_END)/sizeof(struct locations);
	fprintf(stderr,"station,event,phase= %d %d %d\n",station,event,phase);
	fprintf(stderr,"nstat,nevn,npha = %d %d %d\n", nstat, nevn, npha);
	fprintf(stderr,"xmin,xmax,ymin,ymax=%f %f %f %f\n",xmin,xmax,ymin,ymax);

//========================================see how many rays===============
/*
xmesh=16;
ymesh=16;
nx=(xmax-xmin)*xmesh;
ny=(ymax-ymin)*ymesh;
FILE *fp3;
fp3=fopen("../raysnumber/rays number16.txt","w2");
for(i=0;i<npha;i++)
	{
	 if(phase[i].quality<=0.0) continue;
	 iin++;		 	
	 traceit(phase[i].slat-ymin, phase[i].slon-xmin, phase[i].rlat-ymin, phase[i].rlon-xmin, 
			phase[i].scale,xpoint,ypoint,xcell,ycell,lencell,&ncell,&length);

		
		for(j=0;j<ncell;j++) 
		 {
		 jj=xcell[j]+ycell[j]*nx;
		 number[jj]++;	
		 }
	}
fprintf(stderr,"%d %d\n",nx,ny);
for(jj=0;jj<nx*ny;jj++)
{
if (number[jj]>150)
fprintf(fp3,"jj number[jj]=%d %d\n",jj,number[jj]);
//if (number[jj]>1000)
//fprintf(fp3,"------------------------------------------------------------\n");
}
fclose(fp3);
*/
//======================================================test==========================================

grid_time=1;
ifany=1;
FILE *fp;
fp=fopen("../output/raysnum.txt","w2");
FILE *fp17;
fp17=fopen("../output/raysnum2.txt","w2");
//=====================================decide grid time===============================
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
	for(i=0;i<npha;i++)
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
/*	if(xmesh==2)
		{
		  for(jj=0;jj<nx*ny;jj++)
	
	  fprintf(fp17,"%d %d %d\n",number[jj],nx,ny);
		}
*/
	for(jj=0;jj<nx*ny&&aa==0;jj++)
	 {
	 	ifany=0;
	 	if (number[jj]>=numax)
	 	{
		 fprintf(stderr,"jj number[jj] %d %d %f %d\n",jj,number[jj],xmesh,nx);//see the first number >=numax
		 aa=1;
		 ifany=1;		
		 grid_time++;		
		continue;
		}
	}	
}
fclose(fp17);
if(grid_time==1)
{
fprintf(stderr,"raysnumber of every_grid<numx\n");
return;
}
fprintf(stderr,"grid_time=%d \n",grid_time);
array_num=(xmax-xmin)*pow(2,grid_time*2)*(ymax-ymin)*xmesh0*ymesh0/4.0;
//int raysnum[grid_time-1][num];
raysnum=(int**)malloc(sizeof(int*)*(grid_time));
for(i=0;i<grid_time;i++)
raysnum[i]=(int*)malloc(sizeof(int)*array_num);
fprintf(stderr,"arraynum=%d\n",array_num);
//===================================================================================
	nx=(xmax-xmin)*pow(2,grid_time)*xmesh0/2.0;
	ny=(ymax-ymin)*pow(2,grid_time)*ymesh0/2.0;
fprintf(stderr,"nx %d ny %d\n",nx,ny);
	/* allocate storage spaces */
	sol = (struct solution *) malloc(nx*ny * sizeof(struct solution));
	phase2 = (struct rec *) malloc(lseek(fdrec,0L,SEEK_END));
	fprintf(stderr,"sol= %d\n", sol);
	x = (double **) malloc((3*nx*ny+nstat+nevn)*8);
	v = (double **) malloc((3*nx*ny+nstat+nevn)*8);
	u = (double **) malloc((npha+3*nx*ny+1)*8);
	fprintf(stderr,"x,v,u= %d %d %d\n",x,v,u);
	xtemp = (double *) malloc((3*nx*ny+nstat+nevn)*8);
	vtemp = (double *) malloc((3*nx*ny+nstat+nevn)*8);
	utemp = (double *) malloc((npha+3*nx*ny+1)*8);
	fprintf(stderr,"xtemp,vtemp,utemp= %d %d %d\n",xtemp,vtemp,utemp);
	
	/* zero everything here */
	for(i=0;i<nstat;i++) station[i].delay=station[i].wei=station[i].sum=station[i].sum2=0.0;
	for(i=0;i<nevn;i++) event[i].delay=event[i].wei=event[i].sum=event[i].sum2=0.0;
	for(i=0;i<nx*ny;i++) {
		sol[i].ds=sol[i].cos=sol[i].sin=0.0;
		sol[i].wei=sol[i].weicos=sol[i].weisin=0.0;
		sol[i].sum=sol[i].sumcos=sol[i].sumsin=0.0;
		sol[i].sum2=sol[i].sum2cos=sol[i].sum2sin=sol[i].sum2cossin=0.0;
		sol[i].nu=0;
	}
	for(i=0;i<nstat+nevn+3*nx*ny;i++) xtemp[i]=vtemp[i]=0.0;
	for(i=0;i<npha+3*nx*ny;i++) utemp[i]=0.0;

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
		
	 traceit(phase[i].slat-ymin, phase[i].slon-xmin, phase[i].rlat-ymin, phase[i].rlon-xmin, 		phase[i].scale,xpoint,ypoint,xcell,ycell,lencell,&ncell,&length);		
		for(j=0;j<ncell;j++) 
		 {
		 jj=xcell[j]+ycell[j]*nx;
		 number0[jj]++;			
		 raysnum[i2][jj]=number0[jj];
		 }
	}
	
}

//==========================================output========================================
for(i2=0;i2<grid_time;i2++)
{
	fprintf(fp,"---------------------------------------------i2=%d------------------------------------------\n",i2);
	for(jj=0;jj<nx*ny;jj++)	
	{
	if(raysnum[i2][jj]>0)
	fprintf(fp,"jj rays %d %d\n",jj,raysnum[i2][jj]);
	
	}
}
//==========================================encrypt===================================
FILE *fp4;
fp4=fopen("../output/trace.txt","w2");
fprintf(fp4,"j3 xcell ycell lencell   grid_num\n");
FILE *fp5;
fp5=fopen("../output/grid1.txt","w2");	
fprintf(fp5,"j  xcell ycell  lencell   raysnum\n");
FILE *fp7;
fp7=fopen("../output/encrypt.txt","w2");
FILE *fp9;
fp9=fopen("../output/check.txt","w2");
//========initial value=================
for(i=0;i<npha;i++)
{	
	if(phase[i].quality<=0.0) continue;
	xmesh=2*xmesh0/2.0;
	ymesh=2*ymesh0/2.0;
	nx=(xmax-xmin)*xmesh;
	ny=(ymax-ymin)*ymesh;
	FILE *fp6;
	fp6=fopen("../output/poor_quality.txt","w2");	
	if(phase[i].quality<=0.0) 
	{
	fprintf(fp6,"quality <0 rays %d\n",i);
	continue;
	}
	fclose(fp6);
	
	traceit(phase[i].slat-ymin, phase[i].slon-xmin, phase[i].rlat-ymin, phase[i].rlon-xmin,phase[i].scale,xpoint0,ypoint0,xcell0,ycell0,lencell0,&ncell0,&length);
		//========start encrypt=================
	
//======output test==============================================================================================================
			
		fprintf(fp5,"------------------------------------------grid 2*2 for ray %d ncell=%d\n",i,ncell0);
				
		for(j=0;j<ncell0;j++)
		{
		
		jj=xcell0[j]+ycell0[j]*nx;
				
		fprintf(fp5,"%d   %d    %d     %f   %d\n",j,xcell0[j],ycell0[j],lencell0[j],raysnum[0][jj]);
		}	
		
//===============================================================================================================================		
		for(i2=0;i2<grid_time-1;i2++)
		{
		  xmesh=(float)pow(2,i2+1)*xmesh0/2.0;
		  ymesh=(float)pow(2,i2+1)*ymesh0/2.0;
		  nx=(xmax-xmin)*xmesh;
		  ny=(ymax-ymin)*ymesh;
		  ncell_sdr=0;
		  for(j=0;j<ncell0;j++)
		  {
		  jj=xcell0[j]+ycell0[j]*nx;	
			if(raysnum[i2][jj]>=numax)
			{
			fprintf(fp7,"\n-------------ray %d for encrypt %d time------ncell=%d----------\n",i,i2+1,ncell0);			
			xx1=xpoint0[j];yy1=ypoint0[j];
			xx2=xpoint0[j+1];yy2=ypoint0[j+1];					
			traceagain(yy1,xx1,yy2,xx2,phase[i].scale,xpoint_temp,ypoint_temp,xcell_temp,ycell_temp,lencell_temp,&nncell,&length);
			fprintf(fp7,"j=%d jj=%d raysnumbefore=%d encrypt_num=%d \n",j,jj,raysnum[i2][jj],nncell);
				for(j2=0;j2<nncell;j2++)
 				{
				if(lencell_temp[j2]>0.0000001)
				{
				xpoint_sdr[ncell_sdr]=xpoint_temp[j2];
				ypoint_sdr[ncell_sdr]=ypoint_temp[j2];
				xcell_sdr[ncell_sdr]=xcell_temp[j2];
				ycell_sdr[ncell_sdr]=ycell_temp[j2];
				lencell_sdr[ncell_sdr]=lencell_temp[j2];			
				jjnew=xcell_sdr[ncell_sdr]+ycell_sdr[ncell_sdr]*nx*2;				
				fprintf(fp7,"-------- jjnew=%d j=%d %d %d %f raysnum=%d \n",jjnew,ncell_sdr,xcell_sdr[ncell_sdr],ycell_sdr[ncell_sdr],lencell_sdr[ncell_sdr],raysnum[i2+1][jjnew]);
				ncell_sdr++;
				
				}
				}
				xpoint_sdr[ncell_sdr]=xpoint_temp[nncell];
				ypoint_sdr[ncell_sdr]=ypoint_temp[nncell];
			}
			else
			{
			xpoint_sdr[ncell_sdr]=xpoint0[j];
			ypoint_sdr[ncell_sdr]=ypoint0[j];
			xcell_sdr[ncell_sdr]=xcell0[j]*2;
			ycell_sdr[ncell_sdr]=ycell0[j]*2;
			lencell_sdr[ncell_sdr]=lencell0[j];
			jjnew=xcell_sdr[ncell_sdr]+ycell_sdr[ncell_sdr]*nx*2;	
ncell_sdr++;
			xpoint_sdr[ncell_sdr]=xpoint0[j+1];
			ypoint_sdr[ncell_sdr]=ypoint0[j+1];
			}
		  }
		  //=============transfer data==============
		  for(j3=0;j3<ncell_sdr;j3++)
		  {
		  xpoint0[j3]=xpoint_sdr[j3];
		  ypoint0[j3]=ypoint_sdr[j3];
		  xcell0[j3]=xcell_sdr[j3];
		  ycell0[j3]=ycell_sdr[j3];
		  lencell0[j3]=lencell_sdr[j3];
		  }
		  xpoint0[ncell_sdr]=xpoint_sdr[ncell_sdr];
		  ypoint0[ncell_sdr]=ypoint_sdr[ncell_sdr];
		  ncell0=ncell_sdr;
		}
		//==========end encrpt=====================
	
	fprintf(fp4,"---------------------------------ray %d ncell_sdr=%d\n",i,ncell_sdr);
	fprintf(fp9,"---------------------------------------%d\n",i);
	
	for(j3=0;j3<ncell_sdr;j3++)
	{
	
	jj=xcell0[j3]+ycell0[j3]*nx*2;
	
	fprintf(fp9,"%d	\n",jj);
	fprintf(fp4,"%d   %d   %d   %f  %d\n",j3,xcell0[j3],ycell0[j3],lencell0[j3],jj);
	sol[jj].wei += lencell0[j3]*lencell0[j3];
	sol[jj].weicos += lencell0[j3]*cos(2.0*PIO180*phase[i].azse)*lencell0[j3]*cos(2.0*PIO180*phase[i].azse);
	sol[jj].weisin += lencell0[j3]*sin(2.0*PIO180*phase[i].azse)*lencell0[j3]*sin(2.0*PIO180*phase[i].azse);

	sol[jj].nu++;
	}
	
}
fclose(fp9);
fclose(fp7);
fclose(fp5);
fclose(fp4);
fclose(fp);
xmesh=pow(2,grid_time)*xmesh0/2.0;
ymesh=pow(2,grid_time)*ymesh0/2.0;
nx=(xmax-xmin)*xmesh;
ny=(ymax-ymin)*ymesh;
//====================================================================================================
	/*get preconditioners for the conjugate gradient*/
	/*assume you already know station & event counts*/
	
	/*get estimate of rank via trace(ata/diag(ata+dtd) = anorm*/
	/*the 20 is from 4**2+1+1+1+1 for the lapacian damping*/
	anorm=0.0;
	for(i=0;i<nstat;i++) if(station[i].count!=0) anorm+=1;
	for(i=0;i<nevn;i++) if(event[i].count!=0) anorm+=1;
FILE *fp8;
fp8=fopen("../output/solnu.txt","w2");
/*
	for(i=0;i<nx*ny;i++) 
	{
	if(sol[i].nu!=0&&sol[i].nu>3000) 
	{
	fprintf(fp8,"%d %d\n",i,sol[i].nu);
	anorm += sol[i].wei/(sol[i].wei+weight*weight*20);
	}
	}
*/
	for(i=0;i<nx*ny;i++) 
	{
	if(sol[i].nu>0) 
	{
	if(raysnum[0][(i%nx)/2+(nx/2)*((i/nx)/2)]>=numax)
	sol[i].nu=raysnum[1][i];
	
	fprintf(fp8,"%d %d\n",i,sol[i].nu);
	anorm += sol[i].wei/(sol[i].wei+weight*weight*20);
	}
	}
fclose(fp8);
#if ANIWEI
	for(i=0;i<nx*ny;i++) if(sol[i].nu!=0) {
		anorm += sol[i].weicos/(sol[i].weicos+aniwei*aniwei*20);
		anorm += sol[i].weisin/(sol[i].weisin+aniwei*aniwei*20);
	}
#endif
	fprintf(stderr,"\nctomo: anorm= %f\n\n",anorm);
	
	
#if ANIWEI
	for(i=0;i<nx*ny;i++) sol[i].weicos = sqrt(sol[i].weicos+20*aniwei*aniwei);
	for(i=0;i<nx*ny;i++) sol[i].weisin = sqrt(sol[i].weisin+20*aniwei*aniwei);
#endif
	for(i=0;i<nx*ny;i++) sol[i].wei = sqrt(sol[i].wei+20*weight*weight);
	for(i=0;i<nstat;i++) station[i].wei = sqrt((float)station[i].count);
	for(i=0;i<nevn;i++) event[i].wei = sqrt((float)event[i].count);
	fprintf(stderr,"weights found\n");
	
	/*figure out how many parameters*/
	for(i=0,ii=0;i<nstat;i++) if(station[i].count!=0) ii++;
	nstat2=ii;
	for(i=0,ii=0;i<nevn;i++) if(event[i].count!=0) ii++;
	nevn2=ii;
	for(i=0,ii=0;i<nx*ny;i++) if(sol[i].nu!=0) ii++;
	nsol2=ii;
	fprintf(stderr,"total stations,events,sols = %d %d %d\n",nstat2,nevn2,nsol2);
	
	/*assign spots to x pointer*/
	
	/*stations*/
	for(i=0,ii=0;i<nstat;i++) {
		if(station[i].count==0) continue;
		x[ii]= &xtemp[i];
		v[ii]= &vtemp[i];
		ii++;
	}
fprintf(stderr,"ii1=%d\n",ii);
	/*event*/
	for(i=0,ii=ii;i<nevn;i++) {
		if(event[i].count==0) continue;
		x[ii]= &xtemp[i+nstat];
		v[ii]= &vtemp[i+nstat];
		
		ii++;
	}
fprintf(stderr,"ii2=%d\n",ii);
	/*slowness*/
	for(i=0,ii=ii;i<nx*ny;i++) {
		if(sol[i].nu==0) continue;
		x[ii]= &xtemp[i+nstat+nevn];
		v[ii]= &vtemp[i+nstat+nevn];		
		ii++;
	}
fprintf(stderr,"ii3=%d\n",ii);
#if ANIWEI
	
	/*cos*/
	for(i=0,ii=ii;i<nx*ny;i++) {
		if(sol[i].nu==0) continue;
		
		x[ii]=&xtemp[i+nstat+nevn+nx*ny];
		v[ii]=&vtemp[i+nstat+nevn+nx*ny];
		ii++;
	}
fprintf(stderr,"ii4=%d %d %d\n",ii,nx,ny);
	/*sin*/
	for(i=0,ii=ii;i<nx*ny;i++) {
		if(sol[i].nu==0) continue;
		x[ii]=&xtemp[i+nstat+nevn+2*nx*ny];
		v[ii]=&vtemp[i+nstat+nevn+2*nx*ny];
		ii++;
	}
fprintf(stderr,"ii5=%d\n",ii);
#endif	
	nparams=ii;
	
	/*put dtimes in u pointer*/
	for(i=0,ii=0;i<npha;i++) {
		if(phase[i].quality<=0.0) continue;
		u[ii] = &utemp[i];
		utemp[i] = phase[i].dtime;
		
		ii++;
	}
	ntimes=ii;
	fprintf(stderr,"number of good times = %d\n",ii);
	/*put in slowness damp eqns*/
	for(ii=ii,i=0;i<nx*ny;i++) {
		if(sol[i].nu==0) continue;
		u[ii] = &utemp[i+npha];
		utemp[i+npha] = 0.0;
		ii++;
	}
#if ANIWEI
	
	/*put in cos damp eqns*/
	for(ii=ii,i=0;i<nx*ny;i++) {
		if(sol[i].nu==0) continue;
		u[ii] = &utemp[i+npha+nx*ny];
		utemp[i+npha+nx*ny] = 0.0;
		ii++;
	}	
	/*put in sin damp eqns*/
	for(ii=ii,i=0;i<nx*ny;i++) {
		if(sol[i].nu==0) continue;
		u[ii] = &utemp[i+npha+2*nx*ny];
		
		utemp[i+npha+2*nx*ny] = 0.0;
		ii++;
	}
#endif
	
	/* now add an eqn to zero station delays*/
#if ANIWEI
	u[ii] = &utemp[1+npha+3*nx*ny];
	utemp[1+npha+3*nx*ny] = 0.0;
#else
	u[ii] = &utemp[1+npha+nx*ny];
	utemp[1+npha+nx*ny] = 0.0;
#endif	
	ii++;
	
	ntimes=ii;
	
	/*now do LSQR*/
	fprintf(stderr,"starting first lsqr,ntimes,nparams= %d %d \n",ntimes,nparams);
	/******************************************************/
	lsqr(ntimes, nparams, x, u, v, nlsqr);
	//===================================add
	/******************************************************/
/*
for(i=0;i<nx*ny;i++)  if(sol[i].nu!=0) 
{
if(xtemp[i+nstat+nevn+2*nx*ny]==0)
fprintf(stderr,"add3===============xtemp=0\n");
}
*/	
	/*undo the preconditioning*/
	for(i=0;i<nstat;i++) if(station[i].count!=0) station[i].delay = xtemp[i]/station[i].wei;
	for(i=0;i<nevn;i++) if(event[i].count!=0) event[i].delay = xtemp[i+nstat]/event[i].wei;
	for(i=0;i<nx*ny;i++) if(sol[i].nu!=0) 	sol[i].ds = xtemp[i+nstat+nevn]/sol[i].wei;
#if ANIWEI
	FILE *fp16;
	fp16=fopen("../output/xtemp0.txt","w2");
fprintf(stderr,"==========================================0705nx=%d ny=%d\n",nx,ny);
	for(i=0;i<nx*ny;i++) if(sol[i].nu!=0) 
	{
	sol[i].cos = xtemp[i+nstat+nevn+nx*ny]/sol[i].weicos;
	if ( xtemp[i+nstat+nevn+nx*ny]==0)
	fprintf(fp16,"%d %d %d %d %d %d %d\n",i,i+nstat+nevn+nx*ny,i,nstat,nevn,nx,ny);
	if ( xtemp[i+nstat+nevn+nx*ny]==0)
	fprintf(fp16,"-----------------------------------\n");
	}
fclose(fp16);
//-----------------------------------------------------------------------------------!=0
for(i=0;i<nx*ny;i++) if(sol[i].nu!=0) 
sol[i].sin = xtemp[i+nstat+nevn+2*nx*ny]/sol[i].weisin;

#endif
	


/*
for(i=0;i<nx*ny;i++) if(sol[i].nu!=0) 
{
//fprintf(stderr,"%d %f\n"i,sol[i].sin)
if(sol[i].sin==0)
fprintf(stderr,"add1===============sol.son=0\n");
}
*/
	
	/* update travel time for slowness & compute residuals */
	fprintf(stderr,"start residual computations\n");
	iin=res1=res2=res3=res4 = 0.0;
	FILE *fp40;
	fp40=fopen("../output/res.txt","w2");
//=====================================================================traceit2=============================
for(i=0,iin=0;i<npha;i++)
{	
	if(phase[i].quality<=0.0) continue;
	xmesh=2*xmesh0/2.0;
	ymesh=2*ymesh0/2.0;
	nx=(xmax-xmin)*xmesh;
	ny=(ymax-ymin)*ymesh;
	
	iin++;		
	traceit(phase[i].slat-ymin, phase[i].slon-xmin, phase[i].rlat-ymin, phase[i].rlon-xmin,phase[i].scale,xpoint0,ypoint0,xcell0,ycell0,lencell0,&ncell0,&length);
		//========start encrypt=================
		for(i2=0;i2<grid_time-1;i2++)
		{
		  xmesh=pow(2,i2+1)*xmesh0/2.0;
		  ymesh=pow(2,i2+1)*ymesh0/2.0;
		  nx=(xmax-xmin)*xmesh;
		  ny=(ymax-ymin)*ymesh;
		  ncell_sdr=0;
		  for(j=0;j<ncell0;j++)
		  {
		  jj=xcell0[j]+ycell0[j]*nx;  
  		 	if(raysnum[i2][jj]>=numax)
			{	
			xx1=xpoint0[j];yy1=ypoint0[j];
			xx2=xpoint0[j+1];yy2=ypoint0[j+1];			
				
			traceagain(yy1,xx1,yy2,xx2,phase[i].scale,xpoint_temp,ypoint_temp,xcell_temp,ycell_temp,lencell_temp,&nncell,&length);
			for(j2=0;j2<nncell;j2++)
 			{
				if(lencell_temp[j2]>0.0000001)
				{
				xpoint_sdr[ncell_sdr]=xpoint_temp[j2];
				ypoint_sdr[ncell_sdr]=ypoint_temp[j2];
				xcell_sdr[ncell_sdr]=xcell_temp[j2];
				ycell_sdr[ncell_sdr]=ycell_temp[j2];
				lencell_sdr[ncell_sdr]=lencell_temp[j2];
				ncell_sdr++;
				}
				}
				xpoint_sdr[ncell_sdr]=xpoint_temp[nncell];
				ypoint_sdr[ncell_sdr]=ypoint_temp[nncell];

			}
			else
			{
			xpoint_sdr[ncell_sdr]=xpoint0[j];
			ypoint_sdr[ncell_sdr]=ypoint0[j];
			xcell_sdr[ncell_sdr]=xcell0[j]*2;
			ycell_sdr[ncell_sdr]=ycell0[j]*2;
			lencell_sdr[ncell_sdr]=lencell0[j];
			jjnew=xcell_sdr[ncell_sdr]+ycell_sdr[ncell_sdr]*nx*2;
			ncell_sdr++;
			xpoint_sdr[ncell_sdr]=xpoint0[j+1];
			ypoint_sdr[ncell_sdr]=ypoint0[j+1];
			}
		  }
		  //=============transfer data==============
		  for(j3=0;j3<ncell_sdr;j3++)
		  {
		  xpoint0[j3]=xpoint_sdr[j3];
		  ypoint0[j3]=ypoint_sdr[j3];
		  xcell0[j3]=xcell_sdr[j3];
		  ycell0[j3]=ycell_sdr[j3];
		  lencell0[j3]=lencell_sdr[j3];
		  }
		  xpoint0[ncell_sdr]=xpoint_sdr[ncell_sdr];
		  ypoint0[ncell_sdr]=ypoint_sdr[ncell_sdr];
		  ncell0=ncell_sdr;
		}
		//==========end encrpt=====================
	time=phase[i].dtime;
	fprintf(fp40,"%f \n",time);
	for(j3=0;j3<ncell_sdr;j3++)
	{
	jj=xcell0[j3]+ycell0[j3]*nx*2;
	time -= lencell0[j3]*sol[jj].ds;
	fprintf(fp40,"%f ",time);
#if ANIWEI
			time -= lencell0[j3]*sol[jj].cos*cos(2.0*PIO180*phase[i].azse);
			time -= lencell0[j3]*sol[jj].sin*sin(2.0*PIO180*phase[i].azse);
#endif
			
	}
//fprintf(fp40,"%f %f %f %f\n",time,station[phase[i].stnno].delay,event[phase[i].evnno].delay,phase[i].quality);
	time -= station[phase[i].stnno].delay;			
	time -= event[phase[i].evnno].delay;
	fprintf(fp40,"%f \n",time);			
	phase[i].dtime=time;
	wei = phase[i].quality;
	res1 += fabs(time) * wei;
	res2 += time * time * wei;
		
	res3 += time * wei;
	res4 += wei;
	
}
fclose(fp40);
xmesh=pow(2,grid_time)*xmesh0/2.0;
ymesh=pow(2,grid_time)*ymesh0/2.0;
nx=(xmax-xmin)*xmesh;
ny=(ymax-ymin)*ymesh;
/*
for(i=0;i<nx*ny;i++) if(sol[i].nu!=0) 
{

if(sol[i].sin==0)
fprintf(stderr,"add2===============sol.son=0 nx ny%d %d %d\n",nx,ny,i);
}*/
//====================================================================================================================	
	fprintf(stderr,"\nres2 w/o damping residuals= %f\n",res2);
	fprintf(stderr,"anorm,log(res2/res4),penalty = %f, %f, %f\n",anorm,log(res2/res4),2*(anorm+1)/(res4-anorm-2));
	fprintf(stderr,"weight,aniwei= %f, %f\n",weight,aniwei);
	fprintf(stderr,"AIC= %f\n", log(res2/res4)+2*(anorm+1)/res4 );
	fprintf(stderr,"AICc= %f\n", log(res2/res4)+2*(anorm+1)/(res4-anorm-2) );
	fprintf(stderr,"FPE= %f\n", (res2/res4)*(res4+anorm)/(res4-anorm) );
	fprintf(stderr,"SIC= %f\n", log(res2/res4)+log(res4)*anorm/res4 );
	fprintf(stderr,"HQ= %f\n", log(res2/res4)+log(log(res4))*anorm/res4 );
	fprintf(stderr,"Cp= %f\n", res2/(.7065)+res4-2*anorm );
	fprintf(stderr,"AICth= %f\n", log(res2/res4)+6*(anorm+1)/res4 );
	fprintf(stderr,"AICs= %f\n\n", log(res2/res4)+(ntimes/nparams)*2*(anorm+1)/res4 );
	
	/*update res2 to account for the damping equations*/
	
	res2slo=res2cos=res2sin=0.0;
	for(j=0;j<ny;j++) for(i=0;i<nx;i++) {
		if(sol[i+nx*j].nu==0) continue;
		/*slo terms*/
		temp = -4.0*weight*sol[i+nx*j].ds;
//--------------------------------------------------------------------------------------
		
		jj2=i/2+(nx/2)*(j/2);
		jj22=(i/2)*2+nx*(j/2)*2;
		
		//=================================================1/2
		if(raysnum[0][jj2]<numax&&raysnum[0][jj2]>=0)
		{
		//left	
		if(i-2>=0 && sol[i+nx*j-1].nu!=0) temp += 0.5*weight*sol[i+nx*j-1].ds + 0.5*weight*sol[i+nx*(j+1)-1].ds;
		if(i-2>=0 && sol[i+nx*j-1].nu==0&& sol[i+nx*j-2].nu!=0) temp += weight*sol[i+nx*j-2].ds;
		//right
		if(i+2<nx && raysnum[0][jj2+1]>=numax) temp += 0.5*weight*sol[i+nx*j+2].ds + 0.5*weight*sol[i+nx*(j+1)+2].ds;
		if(i+2<nx && raysnum[0][jj2+1]>0&&raysnum[0][jj2+1]<numax) temp += weight*sol[i+nx*j+2].ds;
		//down	
		if(j-2>=0 && sol[i+nx*(j-1)].nu!=0) temp += 0.5*weight*sol[i+nx*(j-1)].ds + 0.5*weight*sol[i+nx*(j-1)+1].ds;
		if(j-2>=0 && sol[i+nx*(j-1)].nu==0&& sol[i+nx*(j-2)].nu!=0) temp += weight*sol[i+nx*(j-2)].ds;	
		//up	
		if(j+2<ny && raysnum[0][jj2+nx/2]>=numax) temp += 0.5*weight*sol[i+nx*(j+2)].ds +0.5*weight*sol[i+nx*(j+2)+1].ds;
		if(j+2<ny && raysnum[0][jj2+nx/2]>=0 && raysnum[0][jj2+nx/2]<numax) temp += weight*sol[i+nx*(j+2)].ds;
		}
		//=================================================1/4
		if(raysnum[0][jj2]>=numax)	
		{
		//left		
		if(i-1>=0 && sol[i+nx*j-1].nu!=0) temp += weight*sol[i+nx*j-1].ds;
		if(i-2>=0 && sol[i+nx*j-1].nu==0 && sol[jj22-2].nu!=0) temp += weight*sol[jj22-2].ds;
		//right
		if(i+1<nx && sol[i+nx*j+1].nu!=0) temp += weight*sol[i+nx*j+1].ds;
		if(i+2<nx && sol[i+nx*j+1].nu==0 && sol[jj22+2].nu!=0) temp += weight*sol[jj22+2].ds;
		//down	
		if(j-1>=0 && sol[i+nx*(j-1)].nu!=0) temp += weight*sol[i+nx*(j-1)].ds;
		if(j>=2&&sol[i+nx*(j-1)].nu==0&& sol[jj22-2*nx].nu!=0) temp += weight*sol[jj22-2*nx].ds;			
		//up
		if(j+1<ny && sol[i+nx*(j+1)].nu!=0) temp += weight*sol[i+nx*(j+1)].ds;
		if(j+2<ny && sol[i+nx*(j+1)].nu==0 &&  sol[jj22+2*nx].nu!=0) temp += weight*sol[jj22+2*nx].ds;		
		}
//-------------------------------------------------------------------
		res2slo += temp*temp/(weight*weight);
		res2 += temp*temp;
#if ANIWEI
		/*cos terms*/
		temp = -4.0*aniwei*sol[i+nx*j].cos;
		
//-------------------------------------------------------------------
//=================================================1/2
		if(raysnum[0][jj2]<numax&&raysnum[0][jj2]>=0)
		{
		//left	
		if(i-2>=0 && sol[i+nx*j-1].nu!=0) temp += 0.5*aniwei*sol[i+nx*j-1].cos + 0.5*aniwei*sol[i+nx*(j+1)-1].cos;
		if(i-2>=0 && sol[i+nx*j-1].nu==0&& sol[i+nx*j-2].nu!=0) temp += aniwei*sol[i+nx*j-2].cos;
		//right
		if(i+2<nx && raysnum[0][jj2+1]>=numax) temp += 0.5*aniwei*sol[i+nx*j+2].cos + 0.5*aniwei*sol[i+nx*(j+1)+2].cos;
		if(i+2<nx && raysnum[0][jj2+1]>0&&raysnum[0][jj2+1]<numax) temp += aniwei*sol[i+nx*j+2].cos;
		//down	
		if(j-2>=0 && sol[i+nx*(j-1)].nu!=0) temp += 0.5*aniwei*sol[i+nx*(j-1)].cos + 0.5*aniwei*sol[i+nx*(j-1)+1].cos;
		if(j-2>=0 && sol[i+nx*(j-1)].nu==0&& sol[i+nx*(j-2)].nu!=0) temp += aniwei*sol[i+nx*(j-2)].cos;	
		//up	
		if(j+2<ny && raysnum[0][jj2+nx/2]>=numax) temp += 0.5*aniwei*sol[i+nx*(j+2)].cos +0.5*aniwei*sol[i+nx*(j+2)+1].cos;
		if(j+2<ny && raysnum[0][jj2+nx/2]>=0 && raysnum[0][jj2+nx/2]<numax) temp += aniwei*sol[i+nx*(j+2)].cos;
		}
		//=================================================1/4
		if(raysnum[0][jj2]>=numax)		
		{
		//left		
		if(i-1>=0 && sol[i+nx*j-1].nu!=0) temp += aniwei*sol[i+nx*j-1].cos;
		if(i-2>=0 && sol[i+nx*j-1].nu==0 && sol[jj22-2].nu!=0) temp += aniwei*sol[jj22-2].cos;
		//right
		if(i+1<nx && sol[i+nx*j+1].nu!=0) temp += aniwei*sol[i+nx*j+1].cos;
		if(i+2<nx && sol[i+nx*j+1].nu==0 && sol[jj22+2].nu!=0) temp += aniwei*sol[jj22+2].cos;
		//down	
		if(j-1>=0 && sol[i+nx*(j-1)].nu!=0) temp += aniwei*sol[i+nx*(j-1)].cos;
		if(j>=2&&sol[i+nx*(j-1)].nu==0&& sol[jj22-2*nx].nu!=0) temp += aniwei*sol[jj22-2*nx].cos;			
		//up
		if(j+1<ny && sol[i+nx*(j+1)].nu!=0) temp += aniwei*sol[i+nx*(j+1)].cos;
		if(j+2<ny && sol[i+nx*(j+1)].nu==0 &&  sol[jj22+2*nx].nu!=0) temp += aniwei*sol[jj22+2*nx].cos;		
		}

//-------------------------------------------------------------------------------------------------
		res2cos += temp*temp/(aniwei*aniwei);
		res2 += temp*temp;
		/*sin terms*/
		
		temp = -4.0*aniwei*sol[i+nx*j].sin;
//=================================================1/2
		if(raysnum[0][jj2]<numax&&raysnum[0][jj2]>=0)
		{
		//left	
		if(i-2>=0 && sol[i+nx*j-1].nu!=0) temp += 0.5*aniwei*sol[i+nx*j-1].sin + 0.5*aniwei*sol[i+nx*(j+1)-1].sin;
		if(i-2>=0 && sol[i+nx*j-1].nu==0&& sol[i+nx*j-2].nu!=0) temp += aniwei*sol[i+nx*j-2].sin;
		//right
		if(i+2<nx && raysnum[0][jj2+1]>=numax) temp += 0.5*aniwei*sol[i+nx*j+2].sin + 0.5*aniwei*sol[i+nx*(j+1)+2].sin;
		if(i+2<nx && raysnum[0][jj2+1]>0&&raysnum[0][jj2+1]<numax) temp += aniwei*sol[i+nx*j+2].sin;
		//down	
		if(j-2>=0 && sol[i+nx*(j-1)].nu!=0) temp += 0.5*aniwei*sol[i+nx*(j-1)].sin + 0.5*aniwei*sol[i+nx*(j-1)+1].sin;
		if(j-2>=0 && sol[i+nx*(j-1)].nu==0&& sol[i+nx*(j-2)].nu!=0) temp += aniwei*sol[i+nx*(j-2)].sin;	
		//up	
		if(j+2<ny && raysnum[0][jj2+nx/2]>=numax) temp += 0.5*aniwei*sol[i+nx*(j+2)].sin +0.5*aniwei*sol[i+nx*(j+2)+1].sin;
		if(j+2<ny && raysnum[0][jj2+nx/2]>=0 && raysnum[0][jj2+nx/2]<numax) temp += aniwei*sol[i+nx*(j+2)].sin;
		}
		//=================================================1/4
		if(raysnum[0][jj2]>=numax)		
		{
		//left		
		if(i-1>=0 && sol[i+nx*j-1].nu!=0) temp += aniwei*sol[i+nx*j-1].sin;
		if(i-2>=0 && sol[i+nx*j-1].nu==0 && sol[jj22-2].nu!=0) temp += aniwei*sol[jj22-2].sin;
		//right
		if(i+1<nx && sol[i+nx*j+1].nu!=0) temp += aniwei*sol[i+nx*j+1].sin;
		if(i+2<nx && sol[i+nx*j+1].nu==0 && sol[jj22+2].nu!=0) temp += aniwei*sol[jj22+2].sin;
		//down	
		if(j-1>=0 && sol[i+nx*(j-1)].nu!=0) temp += aniwei*sol[i+nx*(j-1)].sin;
		if(j>=2&&sol[i+nx*(j-1)].nu==0&& sol[jj22-2*nx].nu!=0) temp += aniwei*sol[jj22-2*nx].sin;			
		//up
		if(j+1<ny && sol[i+nx*(j+1)].nu!=0) temp += aniwei*sol[i+nx*(j+1)].sin;
		if(j+2<ny && sol[i+nx*(j+1)].nu==0 &&  sol[jj22+2*nx].nu!=0) temp += aniwei*sol[jj22+2*nx].sin;		
		}

//-------------------------------------------------------------------------------------------------
		res2sin += temp*temp/(aniwei*aniwei);
		res2 += temp*temp;
#endif
	}
	fprintf(stderr,"res2 with damping residuals= %f\n",res2);
	
	fprintf(stderr,"weight,aniwei= %f %f\n",weight,aniwei);
	
	/*print it out*/
	/*res1 & res2 are weighted residuals*/
	fprintf(stderr,"iin,res1,res2, %10d %10.6f %10.6f\n", iin, res1, res2);
	fprintf(stderr,"res3,res4 = %10.6f %10.5f\n", res3, res4);
	fprintf(stderr,"ntimes,nparams= %d, %d\n",ntimes,nparams);
	fprintf(stderr,"rms1 = %10.5f\n",sqrt(res2/(ntimes-nparams)));
	//fprintf(stderr,"rms2 = %10.5f\n",sqrt(res2/iin));
	fprintf(stderr,"average delay= %20.15f\n", (res3 / res4));
	
	
#if !NBOOTV
	/*this uses the lrsq estimates of variance*/
	/*the sin-cos covariances are not estimated*/
	
	/*undo the preconditioning, and assign*/
	for(i=0;i<nstat;i++) if(station[i].count!=0) station[i].sum2 = vtemp[i]/station[i].wei/station[i].wei;
	for(i=0;i<nevn;i++) if(event[i].count!=0) event[i].sum2 = vtemp[i+nstat]/event[i].wei/event[i].wei;
	for(i=0;i<nx*ny;i++) if(sol[i].nu!=0) sol[i].sum2 = vtemp[i+nstat+nevn]/sol[i].wei/sol[i].wei;
#if ANIWEI
	for(i=0;i<nx*ny;i++) if(sol[i].nu!=0) sol[i].sum2cos = vtemp[i+nstat+nevn+nx*ny]/sol[i].weicos/sol[i].weicos;
	for(i=0;i<nx*ny;i++) if(sol[i].nu!=0) sol[i].sum2sin = vtemp[i+nstat+nevn+2*nx*ny]/sol[i].weisin/sol[i].weisin;
#endif
#endif
	
#if NBOOTV
	/*do bootstrap errors, put results into sum & sum2*/
	fprintf(stderr,"\nstart bootstrap resolution\n");
	/*put into sum2 and sin*/
	for(i=0;i<nstat;i++) station[i].sum=station[i].sum2=0.0;
	
	for(i=0;i<nevn;i++) event[i].sum=event[i].sum2=0.0;
	for(i=0;i<nx*ny;i++) sol[i].sum=sol[i].sum2=0.0;
	
	/*store a copy of phase in phase2*/
	for(i=0;i<npha;i++) phase2[i]=phase[i];
	
	/*do the bootstrap iterations*/
	for(k=0;k<nbootv;k++) {
		fprintf(stderr,"bootstrap errors k= %d\n",k);
		
		
		/*fill up the traveltimes in utemp*/
		for(i=0;i<npha;i++) {
			if(phase2[i].quality<=0) continue;
			while(phase2[l=npha*ran1(4321)].quality<=0);
			phase[i] = phase2[l];
			utemp[i]=phase[i].ttime;
		}
		/*put in damp eqns*/
		for(i=0;i<nx*ny;i++) {
			if(sol[i].nu==0) continue;
			utemp[i+npha]=0.0;
#if ANIWEI
			utemp[i+npha+nx*ny]=utemp[i+npha+2*nx*ny]=0.0;
#endif
		}	
		
		/*now do LSQR*/
		fprintf(stderr,"starting second lsqr,ntimes,nparams,u[0]= %d %d %f\n",ntimes,nparams,*u[0]);
		lsqr(ntimes, nparams, x, u, v, nlsqr);
		
		/*undo the preconditioning*/
		for(j=0;j<nstat;j++) if(station[j].count!=0) xtemp[j] = xtemp[j]/station[j].wei;
		for(j=0;j<nevn;j++) if(event[j].count!=0) xtemp[j+nstat] = xtemp[j+nstat]/event[j].wei;
		for(j=0;j<nx*ny;j++) if(sol[j].nu!=0) xtemp[j+nstat+nevn] = xtemp[j+nstat+nevn]/sol[j].wei;
#if ANIWEI
		for(i=0;i<nx*ny;i++) if(sol[i].nu!=0) xtemp[i+nstat+nevn+nx*ny] /= sol[i].weicos;
		for(i=0;i<nx*ny;i++) if(sol[i].nu!=0) xtemp[i+nstat+nevn+2*nx*ny] /= sol[i].weisin;
#endif
		
		/*update the variance in sum*/
		/*stations*/
		for(j=0;j<nstat;j++) if(station[j].count!=0) {
			
			station[j].sum2 += (xtemp[j]*xtemp[j]);
			station[j].sum += xtemp[j];
		}
		/*events*/
		for(j=0;j<nevn;j++) if(event[j].count!=0) {
			event[j].sum2 += (xtemp[j+nstat]*xtemp[j+nstat]);
			event[j].sum += xtemp[j+nstat];
		}
		/*slownesses*/
		for(j=0;j<nx*ny;j++) if(sol[j].nu!=0) {
			sol[j].sum2 += (xtemp[j+nstat+nevn]*xtemp[j+nstat+nevn]);
			sol[j].sum += xtemp[j+nstat+nevn];
		}
#if ANIWEI
		/*cosines*/
		
		for(j=0;j<nx*ny;j++) if(sol[j].nu!=0) {
			sol[j].sum2cos += xtemp[j+nstat+nevn+nx*ny]*xtemp[j+nstat+nevn+nx*ny];
			sol[j].sumcos += xtemp[j+nstat+nevn+nx*ny];
		}
		/*sines*/
		for(j=0;j<nx*ny;j++) if(sol[j].nu!=0) {
			sol[j].sum2sin += xtemp[j+nstat+nevn+2*nx*ny]*xtemp[j+nstat+nevn+2*nx*ny];
			sol[j].sumsin += xtemp[j+nstat+nevn+2*nx*ny];
		}	
		/*cos-sin*/
		for(j=0;j<nx*ny;j++) if(sol[j].nu!=0) {
			sol[j].sum2cossin += xtemp[j+nstat+nevn+nx*ny]*xtemp[j+nstat+nevn+2*nx*ny];
		}	
#endif
	}
for(i=0;i<nx*ny;i++) if(sol[i].nu!=0) 
{
if(sol[i].sin==0)
fprintf(stderr,"add3===============sol.sin=0\n");
}
	/*compute mean in sum & variance in sum2*/
	/*stations*/
	for(j=0;j<nstat;j++) {
		if(station[j].count==0) continue;
		station[j].sum /= nbootv; /*bootstrap estimate ave*/
		station[j].sum2 = (station[j].sum2 - nbootv*station[j].sum*station[j].sum)/(nbootv-1);
	}
	/*events*/
	for(j=0;j<nevn;j++) {
		if(event[j].count==0) continue;
		event[j].sum /= nbootv;
		event[j].sum2 = (event[j].sum2 - nbootv*event[j].sum*event[j].sum)/(nbootv-1);
	}
	/*slowness*/
	for(j=0;j<nx*ny;j++) {
		if(sol[j].nu==0) continue;
		sol[j].sum /= nbootv;
		sol[j].sum2 = (sol[j].sum2 - nbootv*sol[j].sum*sol[j].sum)/(nbootv-1);;
	}
#if ANIWEI
	/*cosines*/
	for(j=0;j<nx*ny;j++) {
		if(sol[j].nu==0) continue;
		sol[j].sumcos /= nbootv;
		sol[j].sum2cos = (sol[j].sum2cos - nbootv*sol[j].sumcos*sol[j].sumcos)/(nbootv-1);;
	}
	/*sines*/
	
	for(j=0;j<nx*ny;j++) {
		if(sol[j].nu==0) continue;
		sol[j].sumsin /= nbootv;
		sol[j].sum2sin = (sol[j].sum2sin - nbootv*sol[j].sumsin*sol[j].sumsin)/(nbootv-1);;
	}
	/*cos-sin*/
	for(j=0;j<nx*ny;j++) {
		if(sol[j].nu==0) continue;
		sol[j].sum2cossin = (sol[j].sum2cossin - nbootv*sol[j].sumcos*sol[j].sumsin)/(nbootv-1);;
	}
#endif
	fprintf(stderr,"sum of variances = %g\n",time);
#endif
	/*end of NBOOTV if*/
	/* write it all out now */
	fprintf(stderr,"start writing it out\n");
	i = write(fdsol, sol, nx * ny * sizeof(struct solution));
	fprintf(stderr,"sol write= %d\n", i);	
	/* all done */
	
	fprintf(stderr,"ctomo all done now\n\n");
	FILE *fp15;
	fp15=fopen("../output/solb.txt","w2");
	for(jj=0;jj<nx*ny;jj++)
	{
	if(sol[jj].nu==0) continue;
	fprintf(fp15,"%f %f %f %f %f %f %f %f %f %f %f %f %f %d\n",sol[jj].ds,sol[jj].cos,sol[jj].sin,sol[jj].wei,sol[jj].weicos,sol[jj].weisin,sol[jj].sum,sol[jj].sumcos,sol[jj].sumsin,sol[jj].sum2,sol[jj].sum2cos,sol[jj].sum2sin,sol[jj].sum2cossin,sol[jj].nu);
	if(sol[jj].sin==0&&sol[jj].cos==0)
	fprintf(fp15,"---------------------------------------\n");
}
	for(i=0;i<grid_time;i++)
	free(raysnum[i]);
	free(raysnum);
}

/*******************************************************************/
/*****************************************************************/
/*these are the traceing subroutines*/

//=====================================================================\D2\D4\CF\C2Ϊ\D0޸Ĳ\BF\B7\D6================================================================================/
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

/*********************traceagain***************************************/
traceagain(slat,slon, rlat, rlon, scale,xpoint,ypoint, xcell,ycell,lencell,nncell,length)
/* computes the travel time from the locations */
/* assumes straight lines */
/*xcell,ycell,are matrixes that contain the cell number*/
/**ncell is the total number of cells the ray crosses*/
float	slat, slon, rlat, rlon, scale;
int *nncell;
int *xcell,*ycell;
double *lencell,*xpoint,*ypoint,*length;
{	
	double	x1, y1, x2, y2, dx, dy;
	double	len, x, y;
	int		ix, iy;
	double	xnew, ynew;
	int		ixnew, iynew;
	double	alpha, sgnalp;
	double	amax();
	/*zero some stuff*/
	*nncell=0;
	
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
	mesh=xmesh*2;	
	nnx=(xmax-xmin)*mesh;
	nny=(ymax-ymin)*mesh;
	/* set some parameters */
	dx = x2 - x1;
	dy = y2 - y1;
	if(dy == 0.0) dy = 0.0001;
	*length = sqrt(dx*dx + dy*dy) * scale;
	ix = intt(x1 * mesh);
	iy = intt(y1 * mesh);
	alpha = dx / dy;
	sgnalp = sign(alpha);
	
	/* now trace */
	while (x<x2) {
		ynew = (iy + amax(sgnalp, 0.0)) / mesh;
		xnew = x + alpha * (ynew - y);
		if(xnew>x2) {
			xnew = x2;
			ynew = y2;
		}
		ixnew = intt(xnew * mesh);
		iynew = iy + sgnalp;
		if(ixnew>ix) {
			ixnew = ix + 1;
			xnew = ((float) ixnew) / mesh;
			ynew = y + (xnew - x) / alpha;
			iynew = iy;
		}
		dx = xnew - x;
		dy = ynew - y;		
		len = sqrt(dx*dx + dy*dy) * scale;
		
		/* this is the important loop */
		
		if(ix >= 0 && ix<nnx && iy >= 0 && iy<nny) {
			xcell[*nncell]=ix;
			ycell[*nncell]=iy;
			xpoint[*nncell]=x;
			ypoint[*nncell]=y;
			lencell[*nncell]=len;			
			(*nncell)++;			
		}
		ix = ixnew;
		iy = iynew;
		x = xnew;
		y = ynew;
		xpoint[*nncell]=x;
		ypoint[*nncell]=y;		
		
	}
	return;
}
/*********************lsqr routines *********************************/
lsqr(m, n, x, u, v,itmax)
/* subroutine to solve the linear tomographic problem Ax=u using */
/* the lsqr algorithm from G.Nolet, Seismic Tomography, p18 */
/*this version has more error control stuff in it*/
/* m is number of data */
/* n is number of unknown */
/* x(n) is the solution */
/* u(m) is the data (overwritten), passed in as phase*/
/* scratch array v(n) is zeroed and then overwritten by the variances*/
/* the array w(n) & v(n) is dynamically allocated */
/* avpu(m,n,u,v) computes u=u+A*v, ie the forward problem */
/* atupv(m,n,u,v) computes v=v+At*u, ie the backprojection*/
int			m, n;
int			itmax;
double		**x,**u,**v;
/*note that x and v are pointers to pointers*/
{
	double damp,anorm,acond,rnorm,arnorm,xnorm;
	double	*w;
	double *se;
	
	int		i, j;
	
	double sqrt(), fabs(), normlz();
	
	double alfa,bbnorm,beta;
	double cs,cs1,cs2,dampsq,ddnorm,delta;
	double gamma,gambar,phi,phibar,psi;
	double res1,res2,rho,rhobar,rhbar1,rhbar2,rhs;
	double sn,sn1,sn2,t,tau,test2;
	double theta,t1,t2,t3,xxnorm,z,zbar;
	
	fprintf(stderr,"lsqr: m,n= %d %d\n",m,n);	
	/* define data, soln, & scratch areas*/
	w = (double *) malloc(8*n);
	se = (double *) malloc(8*n);
	
	/*initialize everything that needs it*/
	damp=anorm=acond=bbnorm=ddnorm=res2=xnorm=xxnorm=sn2=z=0.0;
	dampsq=damp*damp;
	cs2= -1.0;
	for(i=0; i<n; i++) *v[i]=*x[i]=se[i]=0.0;
	/*start of lsqr*/
	/*set up first bidiagonalizition vectors*/
	/*these satisfy beta*u=b, alfa*v=atrans*u */
	beta = normlz(m, u);
	atupv(m,n,u,v);
	alfa = normlz(n, v);
	
	for(i=0; i<n; i++) w[i] = *v[i];
	rhobar = alfa;
	phibar = beta;
	rnorm = beta;
	arnorm = alfa*beta;
	test2=alfa/beta;

	fprintf(stderr,"beta,alfa=%f %f\n",beta,alfa);
	fprintf(stderr,"iter,arnorm,rnorm,anorm,acond,xnorm,test2 S %f %f %f %f %f %f\n",
		arnorm,rnorm,anorm,acond,xnorm,test2);
		
	/*the main iteration loop*/
	for (i = 0; i < itmax; i++) {
		/*do the next bidiagonalization for beta,u,alfa,v updates*/
		/*then beta*u=a*v-alfa*u and alfa*v=At*u-beta*v*/
		
		for (j=0;j<m;j++) *u[j] *= -alfa;	/* bidiagonalization */		
		avpu(m,n,u,v);
		beta = normlz(m, u);
		bbnorm += alfa*alfa + beta*beta + dampsq;
		for(j=0;j<n;j++) *v[j] *= -beta;
		atupv(m,n,u,v);
		alfa = normlz(n, v);
		/*use plate rotation to eliminate damping parameter*/
		rhbar2 = rhobar*rhobar + dampsq;
		rhbar1 = sqrt(rhbar2);
		cs1 = rhobar/rhbar1;
		sn1 = damp/rhbar1;
		psi= sn1*phibar;
		phibar *= cs1;
		/*use plane rotation to elimniate the subdiagonal element (beta) */
		rho = sqrt(rhbar2 + beta * beta);	/* modified QR */
		cs = rhbar1 / rho;
		sn = beta / rho;
		theta = sn * alfa;
		rhobar = -cs * alfa;
		phi = cs * phibar;
		phibar = sn * phibar;	/* phibar is the sqrt of the sum-of-squares */
		tau = sn*phi;
		/*update x,w and errors*/
		t1 = phi / rho;
		t2 = -theta / rho;
		t3 = 1.0/rho;
		
	// fprintf(stderr,"add========%f %f \n",t1,t);	
		for (j = 0; j < n; j++) {
			t= w[j];
			*x[j] = t1 * t + (*x[j]);
			
			w[j] = t2 * t + (*v[j]);
			t = t3*t*t3*t;
			//fprintf(stderr,"%f\n",t);
			se[j] += t;
			ddnorm += t;
		}
		/*plane rotate on right to eleminate theta*/
		/*and use to estimate xnorm */
		delta = sn2*rho;
		gambar = -cs2*rho;
		rhs = phi-delta*z;
		zbar = rhs/gambar;
		xnorm = sqrt(xxnorm+zbar*zbar);
		gamma = sqrt(gambar*gambar+theta*theta);
		cs2 = gambar/gamma;
		sn2 = theta/gamma;
		z = rhs/gamma;
		xxnorm += z*z;
		
		/*estimates norms */
		anorm = sqrt(bbnorm);	/*anorm is estimate of frobenius norm*/
		acond = anorm*sqrt(ddnorm);
		res1 = phibar*phibar;
		res2 += psi*psi;
		rnorm = sqrt(res1+res2);	/*estimate of sqrt(norm)*/
		arnorm = alfa*fabs(tau); /*arnorm is Atranspose X, should be 0 if converged*/
		
		/*get test2 from norms (machine limit at convergence)*/
		test2=arnorm/(anorm*rnorm);
//**********************************************************************************************************not come to here
		fprintf(stderr,"iter,arnorm,rnorm,anorm,acond,xnorm,test2 %d %f %f %f %f %f %f\n",
			i,arnorm,rnorm,anorm,acond,xnorm,test2);
		
		/* arnorm is size of  Atrans*r */
		/* rnorm is sqrt of sos of residual vector r */
		/* anorm is norm of A */
		/* acond is condition number */
		/* xnorm is norm of sln vector x */
		/* test2 < machine limit at convergence */
	}
	
	for(i=0;i<n;i++) *v[i] = se[i]*rnorm*rnorm/(m-n);
	
	return;
}

/*************************************************************************/
double normlz(n, x)
int		n;
double	**x;

/* normalizes vector x */
{
	int i,j;
	double ss,s ;
	
	/*fprintf(stderr,"normlz: start, n= %d\n",n);*/
	
	for(i=0,ss=0.0; i<n; i++) ss += (*x[i])*(*x[i]);
	s= sqrt(ss);
	ss=1.0/s;
	for(i=0; i<n; i++) *x[i] *= ss;
	
	/*fprintf(stderr,"normlz: s, s**2= %f %f\n",s,s*s);*/
	return(s);
}

/*********************************************************************/

avpu(m,n,u,v)
int m,n; 
double **u,**v;
/*the forward projection*/
{
	int		i,j,ii,jj,k,kk,l,iii;
	int ncell;
	int xcell[num],ycell[num];
	double lencell[num], azcell[num];
	double length;
	double exp();
	//=============add============
	int ncell0,nncell;
	int number2[num2],number0[num2];
	int xcell0[num],ycell0[num],xcell_temp[num],ycell_temp[num],xcell_sdr[num],ycell_sdr[num];
	double lencell0[num],lencell_temp[num],lencell_sdr[num];
	double xpoint[num],ypoint[num],xpoint0[num],ypoint0[num],xpoint_temp[num],ypoint_temp[num];
	double xpoint_sdr[num],ypoint_sdr[num];
	int jj1,jj2,jj11,jj22;
	//============================
	/*fprintf(stderr,"avpu: start\n");*/
	int **raysnum;
	/*for each ray call trace*/
	//===========================================================================================================

//======================================================test==========================================

grid_time=1;
ifany=1;
//=====================================decide grid time===============================
for(grid_time=1;ifany>0&&grid_time<grid_max;)
{

aa=0;
xmesh=pow(2,grid_time)*xmesh0/2.0;
ymesh=pow(2,grid_time)*ymesh0/2.0;
nx=(xmax-xmin)*xmesh;
ny=(ymax-ymin)*ymesh;
for(jj2=0;jj2<nx*ny;jj2++)
{
number2[jj2]=0;
}
//fprintf(stderr,"mesh,nx,ny=%f %d %d\n-------------\n",xmesh,nx,ny);
    for(i=0;i<npha;i++)
    {
	 if(phase[i].quality<=0.0) continue;	
	 traceit(phase[i].slat-ymin, phase[i].slon-xmin, phase[i].rlat-ymin, phase[i].rlon-xmin, 
			phase[i].scale,xpoint,ypoint,xcell,ycell,lencell,&ncell,&length);
	for(j=0;j<ncell;j++) 
	{
	jj=xcell[j]+ycell[j]*nx;			
	number2[jj]++;
	}
    }
	//fprintf(stderr,"mesh nx ny=%f %d %d\n",xmesh,nx,ny);
	for(jj=0;jj<nx*ny&&aa==0;jj++)
	 {
	 	ifany=0;
	 	if (number2[jj]>=numax)
	 	{
		 aa=1;
		 ifany=1;		
		 grid_time++;	 
		continue;
		}
	}	
}	
//fprintf(stderr,"-----------------------avpu grid_time=%d\n",grid_time);
//int raysnum[grid_time-1][num];
array_num=(xmax-xmin)*pow(2,grid_time*2)*(ymax-ymin)*xmesh0*ymesh0/4.0;
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
//========initial value=================
for(i=0;i<npha;i++)
{
	if(phase[i].quality<=0.0) continue;	
	xmesh=2*xmesh0/2.0;
	ymesh=2*ymesh0/2.0;
	nx=(xmax-xmin)*xmesh;
	ny=(ymax-ymin)*ymesh;		
	traceit(phase[i].slat-ymin, phase[i].slon-xmin, phase[i].rlat-ymin, phase[i].rlon-xmin,phase[i].scale,xpoint0,ypoint0,xcell0,ycell0,lencell0,&ncell0,&length);
		//========start encrypt=================
	for(i2=0;i2<grid_time-1;i2++)
		{
		  xmesh=pow(2,i2+1)*xmesh0/2.0;
		  ymesh=pow(2,i2+1)*ymesh0/2.0;
		  nx=(xmax-xmin)*xmesh;
		  ny=(ymax-ymin)*ymesh;
		  ncell_sdr=0;
		  for(j=0;j<ncell0;j++)
		  {
		  jj=xcell0[j]+ycell0[j]*nx;	
			if(raysnum[i2][jj]>=numax)
			{
			xx1=xpoint0[j];yy1=ypoint0[j];
			xx2=xpoint0[j+1];yy2=ypoint0[j+1];
	traceagain(yy1,xx1,yy2,xx2,phase[i].scale,xpoint_temp,ypoint_temp,xcell_temp,ycell_temp,lencell_temp,&nncell,&length);
			for(j2=0;j2<nncell;j2++)
 				{
				if(lencell_temp[j2]>0.0000001)
				{
				xpoint_sdr[ncell_sdr]=xpoint_temp[j2];
				ypoint_sdr[ncell_sdr]=ypoint_temp[j2];
				xcell_sdr[ncell_sdr]=xcell_temp[j2];
				ycell_sdr[ncell_sdr]=ycell_temp[j2];
				lencell_sdr[ncell_sdr]=lencell_temp[j2];			
				ncell_sdr++;
				}
				}
				xpoint_sdr[ncell_sdr]=xpoint_temp[nncell];
				ypoint_sdr[ncell_sdr]=ypoint_temp[nncell];
			}
			else
			{
			xpoint_sdr[ncell_sdr]=xpoint0[j];
			ypoint_sdr[ncell_sdr]=ypoint0[j];
			xcell_sdr[ncell_sdr]=xcell0[j]*2;
			ycell_sdr[ncell_sdr]=ycell0[j]*2;
			lencell_sdr[ncell_sdr]=lencell0[j];
			ncell_sdr++;
			xpoint_sdr[ncell_sdr]=xpoint0[j+1];
			ypoint_sdr[ncell_sdr]=ypoint0[j+1];
			}
		  }
		  //=============transfer data==============
		  for(j3=0;j3<ncell_sdr;j3++)
		  {
		  xpoint0[j3]=xpoint_sdr[j3];
		  ypoint0[j3]=ypoint_sdr[j3];
		  xcell0[j3]=xcell_sdr[j3];
		  ycell0[j3]=ycell_sdr[j3];
		  lencell0[j3]=lencell_sdr[j3];
		  }
		  xpoint0[ncell_sdr]=xpoint_sdr[ncell_sdr];
		  ypoint0[ncell_sdr]=ypoint_sdr[ncell_sdr];
		  ncell0=ncell_sdr;
		}
		//==========end encrpt=====================
	for(j3=0;j3<ncell_sdr;j3++)
	{
	jj=xcell0[j3]+ycell0[j3]*nx*2;

#if ANIWEI
			utemp[i] += lencell0[j3]*vtemp[jj+nstat+nevn+nx*ny*4]*cos(2.0*PIO180*phase[i].azse)/sol[jj].weicos
				+ lencell0[j3]*vtemp[jj+nstat+nevn+2*nx*ny*4]*sin(2.0*PIO180*phase[i].azse)/sol[jj].weisin
				+ lencell0[j3]*vtemp[jj+nstat+nevn]/sol[jj].wei;
#else
			utemp[i] += lencell0[j3]*vtemp[jj+nstat+nevn]/sol[jj].wei;
#endif
	}
	utemp[i] += vtemp[phase[i].stnno]/station[phase[i].stnno].wei
			+ vtemp[phase[i].evnno+nstat]/event[phase[i].evnno].wei;
}
xmesh=pow(2,grid_time)*xmesh0/2.0;
ymesh=pow(2,grid_time)*ymesh0/2.0;
nx=(xmax-xmin)*xmesh;
ny=(ymax-ymin)*ymesh;
//=============================================================================================
	
	/*zero station delays*/
	for(i=0;i<nstat;i++) {
		if(station[i].count==0) continue;
#if ANIWEI
		utemp[npha+1+3*nx*ny] += vtemp[i]/station[i].wei;
#else	
		utemp[npha+1+nx*ny] += vtemp[i]/station[i].wei;
#endif	
	}
		
	/*damp slowness here*/
	for(j=0;j<ny;j++) for(i=0;i<nx;i++) {
		if(sol[i+nx*j].nu==0) continue;
		utemp[npha+i+nx*j] += -4.0*weight*vtemp[i+nx*j+nstat+nevn]/sol[i+nx*j].wei;
//--------------------------------------------------------------------------------------
		
		jj2=i/2+(nx/2)*(j/2);
		jj22=(i/2)*2+nx*(j/2)*2;
//=================================================1/2
		if(raysnum[0][jj2]<numax&&raysnum[0][jj2]>=0)
		{
		//left	
		if(i-2>=0 && sol[i+nx*j-1].nu!=0) utemp[npha+i+nx*j] += 0.5*weight*vtemp[nstat+nevn+i+nx*j-1]/sol[i+nx*j-1].wei + 0.5*weight*vtemp[nstat+nevn+i+nx*(j+1)-1]/sol[i+nx*(j+1)-1].wei;
		if(i-2>=0 && sol[i+nx*j-1].nu==0&& sol[i+nx*j-2].nu!=0) utemp[npha+i+nx*j] += weight*vtemp[nstat+nevn+i+nx*j-2]/sol[i+nx*j-2].wei;
		//right
		if(i+2<nx && raysnum[0][jj2+1]>=numax) utemp[npha+i+nx*j] += 0.5*weight*vtemp[nstat+nevn+i+nx*j+2]/sol[i+nx*j+2].wei + 0.5*weight*vtemp[nstat+nevn+i+nx*(j+1)+2]/sol[i+nx*(j+1)+2].wei;
		if(i+2<nx && raysnum[0][jj2+1]>0&&raysnum[0][jj2+1]<numax) utemp[npha+i+nx*j] += weight*vtemp[nstat+nevn+i+nx*j+2]/sol[i+nx*j+2].wei;
		//down	
		if(j-2>=0 && sol[i+nx*(j-1)].nu!=0) utemp[npha+i+nx*j] += 0.5*weight*vtemp[nstat+nevn+i+nx*(j-1)]/sol[i+nx*(j-1)].wei + 0.5*weight*vtemp[nstat+nevn+i+nx*(j-1)+1]/sol[i+nx*(j-1)+1].wei;
		if(j-2>=0 && sol[i+nx*(j-1)].nu==0&& sol[i+nx*(j-2)].nu!=0) utemp[npha+i+nx*j] += weight*vtemp[nstat+nevn+i+nx*(j-2)]/sol[i+nx*(j-2)].wei;	
		//up	
		if(j+2<ny && raysnum[0][jj2+nx/2]>=numax) utemp[npha+i+nx*j] += 0.5*weight*vtemp[nstat+nevn+i+nx*(j+2)]/sol[i+nx*(j+2)].wei +0.5*weight*vtemp[nstat+nevn+i+nx*(j+2)+1]/sol[i+nx*(j+2)+1].wei;
		if(j+2<ny && raysnum[0][jj2+nx/2]>=0 && raysnum[0][jj2+nx/2]<numax) utemp[npha+i+nx*j] += weight*vtemp[nstat+nevn+i+nx*(j+2)]/sol[i+nx*(j+2)].wei;
		}
		
//=================================================1/4
		if(raysnum[0][jj2]>=numax)		
		{
		//left		
		if(i-1>=0 && sol[i+nx*j-1].nu!=0) utemp[npha+i+nx*j] += weight*vtemp[i+nx*j-1+nstat+nevn]/sol[i+nx*j-1].wei;
		if(i-2>=0 && sol[i+nx*j-1].nu==0 && sol[jj22-2].nu!=0) utemp[npha+i+nx*j] += weight*vtemp[jj22-2+nstat+nevn]/sol[jj22-2].wei;
		//right
		if(i+1<nx && sol[i+nx*j+1].nu!=0) utemp[npha+i+nx*j] += weight*vtemp[i+nx*j+1+nstat+nevn]/sol[i+nx*j+1].wei;
		if(i+2<nx && sol[i+nx*j+1].nu==0 && sol[jj22+2].nu!=0) utemp[npha+i+nx*j] += weight*vtemp[jj22+2+nstat+nevn]/sol[jj22+2].wei;
		//down	
		if(j-1>=0 )
		{
			if(sol[i+nx*(j-1)].nu!=0) utemp[npha+i+nx*j] += weight*vtemp[i+nx*(j-1)+nstat+nevn]/sol[i+nx*(j-1)].wei;
			if(j>=2) {if ( sol[i+nx*(j-1)].nu==0&& sol[jj22-2*nx].nu!=0) utemp[npha+i+nx*j] += weight*vtemp[jj22-2*nx+nstat+nevn]/sol[jj22-2*nx].wei;}
		}
		//up
		if(j+1<ny && sol[i+nx*(j+1)].nu!=0) utemp[npha+i+nx*j] += weight*vtemp[i+nx*(j+1)+nstat+nevn]/sol[i+nx*(j+1)].wei;
		if(j+2<ny && sol[i+nx*(j+1)].nu==0 &&  sol[jj22+2*nx].nu!=0) utemp[npha+i+nx*j] += weight*vtemp[jj22+2*nx+nstat+nevn]/sol[jj22+2*nx].wei;
		}
	}







#if ANIWEI
	/*damp cos here*/
	for(j=0;j<ny;j++) for(i=0;i<nx;i++) {
		if(sol[i+nx*j].nu==0) continue;
		utemp[npha+i+nx*j+nx*ny] += -4.0*aniwei*vtemp[i+nx*j+nstat+nevn+nx*ny]/sol[i+nx*j].weicos;
//--------------------------------------------------------------------------------------
		
		jj2=i/2+(nx/2)*(j/2);
		jj22=(i/2)*2+nx*(j/2)*2;
		//=================================================1/2
		if(raysnum[0][jj2]<numax&&raysnum[0][jj2]>=0)
		{
		//left	
		if(i-2>=0 && sol[i+nx*j-1].nu!=0) utemp[npha+i+nx*j+nx*ny] += 0.5*aniwei*vtemp[nstat+nevn+nx*ny+i+nx*j-1]/sol[i+nx*j-1].weicos + 0.5*aniwei*vtemp[nstat+nevn+nx*ny+i+nx*(j+1)-1]/sol[i+nx*(j+1)-1].weicos;
		if(i-2>=0 && sol[i+nx*j-1].nu==0&& sol[i+nx*j-2].nu!=0) utemp[npha+i+nx*j+nx*ny] += aniwei*vtemp[nstat+nevn+nx*ny+i+nx*j-2]/sol[i+nx*j-2].weicos;
		//right
		if(i+2<nx && raysnum[0][jj2+1]>=numax) utemp[npha+i+nx*j+nx*ny] += 0.5*aniwei*vtemp[nstat+nevn+nx*ny+i+nx*j+2]/sol[i+nx*j+2].weicos + 0.5*aniwei*vtemp[nstat+nevn+nx*ny+i+nx*(j+1)+2]/sol[i+nx*(j+1)+2].weicos;
		if(i+2<nx && raysnum[0][jj2+1]>0&&raysnum[0][jj2+1]<numax) utemp[npha+i+nx*j+nx*ny] += aniwei*vtemp[nstat+nevn+nx*ny+i+nx*j+2]/sol[i+nx*j+2].weicos;
		//down	
		if(j-2>=0 && sol[i+nx*(j-1)].nu!=0) utemp[npha+i+nx*j+nx*ny] += 0.5*aniwei*vtemp[nstat+nevn+nx*ny+i+nx*(j-1)]/sol[i+nx*(j-1)].weicos + 0.5*aniwei*vtemp[nstat+nevn+nx*ny+i+nx*(j-1)+1]/sol[i+nx*(j-1)+1].weicos;
		if(j-2>=0 && sol[i+nx*(j-1)].nu==0&& sol[i+nx*(j-2)].nu!=0) utemp[npha+i+nx*j+nx*ny] += aniwei*vtemp[nstat+nevn+nx*ny+i+nx*(j-2)]/sol[i+nx*(j-2)].weicos;	
		//up	
		if(j+2<ny && raysnum[0][jj2+nx/2]>=numax) utemp[npha+i+nx*j+nx*ny] += 0.5*aniwei*vtemp[nstat+nevn+nx*ny+i+nx*(j+2)]/sol[i+nx*(j+2)].weicos +0.5*aniwei*vtemp[nstat+nevn+nx*ny+i+nx*(j+2)+1]/sol[i+nx*(j+2)+1].weicos;
		if(j+2<ny && raysnum[0][jj2+nx/2]>=0 && raysnum[0][jj2+nx/2]<numax) utemp[npha+i+nx*j+nx*ny] += aniwei*vtemp[nstat+nevn+nx*ny+i+nx*(j+2)]/sol[i+nx*(j+2)].weicos;
		}
		
//=================================================1/4
		if(raysnum[0][jj2]>=numax)		
		{
		//left		
		if(i-1>=0 && sol[i+nx*j-1].nu!=0) utemp[npha+i+nx*j+nx*ny] += aniwei*vtemp[i+nx*j-1+nstat+nevn+nx*ny]/sol[i+nx*j-1].weicos;
		if(i-2>=0 && sol[i+nx*j-1].nu==0 && sol[jj22-2].nu!=0) utemp[npha+i+nx*j+nx*ny] += aniwei*vtemp[jj22-2+nstat+nevn+nx*ny]/sol[jj22-2].weicos;
		//right
		if(i+1<nx && sol[i+nx*j+1].nu!=0) utemp[npha+i+nx*j+nx*ny] += aniwei*vtemp[i+nx*j+1+nstat+nevn+nx*ny]/sol[i+nx*j+1].weicos;
		if(i+2<nx && sol[i+nx*j+1].nu==0 && sol[jj22+2].nu!=0) utemp[npha+i+nx*j+nx*ny] += aniwei*vtemp[jj22+2+nstat+nevn+nx*ny]/sol[jj22+2].weicos;
		//down	
		if(j-1>=0 )
		{
			if(sol[i+nx*(j-1)].nu!=0) utemp[npha+i+nx*j+nx*ny] += aniwei*vtemp[i+nx*(j-1)+nstat+nevn+nx*ny]/sol[i+nx*(j-1)].weicos;
			if(j>=2) {if ( sol[i+nx*(j-1)].nu==0&& sol[jj22-2*nx].nu!=0) utemp[npha+i+nx*j+nx*ny] += aniwei*vtemp[jj22-2*nx+nstat+nevn+nx*ny]/sol[jj22-2*nx].weicos;}
		}
		//up
		if(j+1<ny && sol[i+nx*(j+1)].nu!=0) utemp[npha+i+nx*j+nx*ny] += aniwei*vtemp[i+nx*(j+1)+nstat+nevn+nx*ny]/sol[i+nx*(j+1)].weicos;
		if(j+2<ny && sol[i+nx*(j+1)].nu==0 &&  sol[jj22+2*nx].nu!=0) utemp[npha+i+nx*j+nx*ny] += aniwei*vtemp[jj22+2*nx+nstat+nevn+nx*ny]/sol[jj22+2*nx].weicos;
		}
	}
	/*damp sin here*/
	for(j=0;j<ny;j++) for(i=0;i<nx;i++) {
		if(sol[i+nx*j].nu==0) continue;
		utemp[npha+i+nx*j+2*nx*ny] += -4.0*aniwei*vtemp[i+nx*j+nstat+nevn+2*nx*ny]/sol[i+nx*j].weisin;
//--------------------------------------------------------------------------------------
		
		jj2=i/2+(nx/2)*(j/2);
		jj22=(i/2)*2+nx*(j/2)*2;
		//=================================================1/2
		if(raysnum[0][jj2]<numax&&raysnum[0][jj2]>=0)
		{
		//left	
		if(i-2>=0 && sol[i+nx*j-1].nu!=0) utemp[npha+i+nx*j+2*nx*ny] += 0.5*aniwei*vtemp[nstat+nevn+2*nx*ny+i+nx*j-1]/sol[i+nx*j-1].weisin + 0.5*aniwei*vtemp[nstat+nevn+2*nx*ny+i+nx*(j+1)-1]/sol[i+nx*(j+1)-1].weisin;
		if(i-2>=0 && sol[i+nx*j-1].nu==0&& sol[i+nx*j-2].nu!=0) utemp[npha+i+nx*j+2*nx*ny] += aniwei*vtemp[nstat+nevn+2*nx*ny+i+nx*j-2]/sol[i+nx*j-2].weisin;
		//right
		if(i+2<nx && raysnum[0][jj2+1]>=numax) utemp[npha+i+nx*j+2*nx*ny] += 0.5*aniwei*vtemp[nstat+nevn+2*nx*ny+i+nx*j+2]/sol[i+nx*j+2].weisin + 0.5*aniwei*vtemp[nstat+nevn+2*nx*ny+i+nx*(j+1)+2]/sol[i+nx*(j+1)+2].weisin;
		if(i+2<nx && raysnum[0][jj2+1]>0&&raysnum[0][jj2+1]<numax) utemp[npha+i+nx*j+2*nx*ny] += aniwei*vtemp[nstat+nevn+2*nx*ny+i+nx*j+2]/sol[i+nx*j+2].weisin;
		//down	
		if(j-2>=0 && sol[i+nx*(j-1)].nu!=0) utemp[npha+i+nx*j+2*nx*ny] += 0.5*aniwei*vtemp[nstat+nevn+2*nx*ny+i+nx*(j-1)]/sol[i+nx*(j-1)].weisin + 0.5*aniwei*vtemp[nstat+nevn+2*nx*ny+i+nx*(j-1)+1]/sol[i+nx*(j-1)+1].weisin;
		if(j-2>=0 && sol[i+nx*(j-1)].nu==0&& sol[i+nx*(j-2)].nu!=0) utemp[npha+i+nx*j+2*nx*ny] += aniwei*vtemp[nstat+nevn+2*nx*ny+i+nx*(j-2)]/sol[i+nx*(j-2)].weisin;	
		//up	
		if(j+2<ny && raysnum[0][jj2+nx/2]>=numax) utemp[npha+i+nx*j+2*nx*ny] += 0.5*aniwei*vtemp[nstat+nevn+2*nx*ny+i+nx*(j+2)]/sol[i+nx*(j+2)].weisin +0.5*aniwei*vtemp[nstat+nevn+2*nx*ny+i+nx*(j+2)+1]/sol[i+nx*(j+2)+1].weisin;
		if(j+2<ny && raysnum[0][jj2+nx/2]>=0 && raysnum[0][jj2+nx/2]<numax) utemp[npha+i+nx*j+2*nx*ny] += aniwei*vtemp[nstat+nevn+2*nx*ny+i+nx*(j+2)]/sol[i+nx*(j+2)].weisin;
		}
//=================================================1/4
		if(raysnum[0][jj2]>=numax)		
		{
		//left		
		if(i-1>=0 && sol[i+nx*j-1].nu!=0) utemp[npha+i+nx*j+2*nx*ny] += aniwei*vtemp[i+nx*j-1+nstat+nevn+2*nx*ny]/sol[i+nx*j-1].weisin;
		if(i-2>=0 && sol[i+nx*j-1].nu==0 && sol[jj22-2].nu!=0) utemp[npha+i+nx*j+2*nx*ny] += aniwei*vtemp[jj22-2+nstat+nevn+2*nx*ny]/sol[jj22-2].weisin;
		//right
		if(i+1<nx && sol[i+nx*j+1].nu!=0) utemp[npha+i+nx*j+2*nx*ny] += aniwei*vtemp[i+nx*j+1+nstat+nevn+2*nx*ny]/sol[i+nx*j+1].weisin;
		if(i+2<nx && sol[i+nx*j+1].nu==0 && sol[jj22+2].nu!=0) utemp[npha+i+nx*j+2*nx*ny] += aniwei*vtemp[jj22+2+nstat+nevn+2*nx*ny]/sol[jj22+2].weisin;
		//down	
		if(j-1>=0 )
		{
			if(sol[i+nx*(j-1)].nu!=0) utemp[npha+i+nx*j+2*nx*ny] += aniwei*vtemp[i+nx*(j-1)+nstat+nevn+2*nx*ny]/sol[i+nx*(j-1)].weisin;
			if(j>=2) {if ( sol[i+nx*(j-1)].nu==0&& sol[jj22-2*nx].nu!=0) utemp[npha+i+nx*j+2*nx*ny] += aniwei*vtemp[jj22-2*nx+nstat+nevn+2*nx*ny]/sol[jj22-2*nx].weisin;}
		}
		//up
		if(j+1<ny && sol[i+nx*(j+1)].nu!=0) utemp[npha+i+nx*j+2*nx*ny] += aniwei*vtemp[i+nx*(j+1)+nstat+nevn+2*nx*ny]/sol[i+nx*(j+1)].weisin;
		if(j+2<ny && sol[i+nx*(j+1)].nu==0 &&  sol[jj22+2*nx].nu!=0) utemp[npha+i+nx*j+2*nx*ny] += aniwei*vtemp[jj22+2*nx+nstat+nevn+2*nx*ny]/sol[jj22+2*nx].weisin;
		}
		
//-------------------------------------------------------------------	
	}
#endif
	for(i=0;i<grid_time;i++)
	free(raysnum[i]);
	free(raysnum);
	/*fprintf(stderr,"avpu done\n");*/
	
}

atupv(m,n,u,v)
int m,n; double **u,**v;
/*backproject celldistance * time*/
/*u is data vector in phase[].empty*/
/*v is temp solution vector is ds2 & delay2*/
/*this is the gradient operator*/
{
	
	int	i,j,ii,jj,k,kk,l,iii;
	int ncell;
	int xcell[num],ycell[num];
	double lencell[num],azcell[num];
	double length;
	double exp();
	//=============add============
	int ncell0,nncell;
	int number[num2],number0[num2];
	int xcell0[num],ycell0[num],xcell_temp[num],ycell_temp[num],xcell_sdr[num],ycell_sdr[num];
	double lencell0[num],lencell_temp[num],lencell_sdr[num];
	double xpoint[num],ypoint[num],xpoint0[num],ypoint0[num],xpoint_temp[num],ypoint_temp[num];
	double xpoint_sdr[num],ypoint_sdr[num];
	int jj1,jj2,jj11,jj22;
	/*fprintf(stderr,"atupv: start\n");*/
	int **raysnum;
	/*backproject the ray*/
	//========================================================================================================
//fprintf(stderr,"add========into atupv\n");
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
		 //fprintf(stderr,"jj number[jj] %d %d\n",jj,number[jj]);//see the first number >=numax
		 aa=1;
		 ifany=1;		
		 grid_time++;	 
		
		 continue;
		 }
		}	
	}
//fprintf(stderr,"add-----------------------atupv grid_time=%d %d\n",grid_time,num);
//int raysnum[grid_time-1][num];
array_num=(xmax-xmin)*pow(2,grid_time*2)*(ymax-ymin)*xmesh0*ymesh0/4.0;
raysnum=(int**)malloc(sizeof(int*)*(grid_time));
for(i=0;i<grid_time;i++)
raysnum[i]=(int*)malloc(sizeof(int)*array_num);
//fprintf(stderr,"add===================88\n");
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
//fprintf(stderr,"add============start gird\n");

//========initial value=================
for(i=0;i<npha;i++)
{
	if(phase[i].quality<=0.0) continue;	
	xmesh=2*xmesh0/2.0;
	ymesh=2*ymesh0/2.0;
	nx=(xmax-xmin)*xmesh;
	ny=(ymax-ymin)*ymesh;
	if(phase[i].quality<=0.0) 
	continue;		
	traceit(phase[i].slat-ymin, phase[i].slon-xmin, phase[i].rlat-ymin, phase[i].rlon-xmin,phase[i].scale,xpoint0,ypoint0,xcell0,ycell0,lencell0,&ncell0,&length);
		//========start encrypt=================
	for(i2=0;i2<grid_time-1;i2++)
		{
		  xmesh=pow(2,i2+1)*xmesh0/2.0;
		  ymesh=pow(2,i2+1)*ymesh0/2.0;
		  nx=(xmax-xmin)*xmesh;
		  ny=(ymax-ymin)*ymesh;
		  ncell_sdr=0;
		  for(j=0;j<ncell0;j++)
		  {
		  jj=xcell0[j]+ycell0[j]*nx;	
			if(raysnum[i2][jj]>=numax)
			{
			xx1=xpoint0[j];yy1=ypoint0[j];
			xx2=xpoint0[j+1];yy2=ypoint0[j+1];
	traceagain(yy1,xx1,yy2,xx2,phase[i].scale,xpoint_temp,ypoint_temp,xcell_temp,ycell_temp,lencell_temp,&nncell,&length);
			for(j2=0;j2<nncell;j2++)
 				{
				{
				xpoint_sdr[ncell_sdr]=xpoint_temp[j2];
				ypoint_sdr[ncell_sdr]=ypoint_temp[j2];
				xcell_sdr[ncell_sdr]=xcell_temp[j2];
				ycell_sdr[ncell_sdr]=ycell_temp[j2];
				lencell_sdr[ncell_sdr]=lencell_temp[j2];			
				ncell_sdr++;
				}
				}
				xpoint_sdr[ncell_sdr]=xpoint_temp[nncell];
				ypoint_sdr[ncell_sdr]=ypoint_temp[nncell];
			}
			else
			{
			xpoint_sdr[ncell_sdr]=xpoint0[j];
			ypoint_sdr[ncell_sdr]=ypoint0[j];
			xcell_sdr[ncell_sdr]=xcell0[j]*2;
			ycell_sdr[ncell_sdr]=ycell0[j]*2;
			lencell_sdr[ncell_sdr]=lencell0[j];
			ncell_sdr++;
			xpoint_sdr[ncell_sdr]=xpoint0[j+1];
			ypoint_sdr[ncell_sdr]=ypoint0[j+1];
			}
		  }
		  //=============transfer data==============
		  for(j3=0;j3<ncell_sdr;j3++)
		  {
		  xpoint0[j3]=xpoint_sdr[j3];
		  ypoint0[j3]=ypoint_sdr[j3];
		  xcell0[j3]=xcell_sdr[j3];
		  ycell0[j3]=ycell_sdr[j3];
		  lencell0[j3]=lencell_sdr[j3];
		  }
		  xpoint0[ncell_sdr]=xpoint_sdr[ncell_sdr];
		  ypoint0[ncell_sdr]=ypoint_sdr[ncell_sdr];
		  ncell0=ncell_sdr;
		}
		//==========end encrpt=====================
	for(j3=0;j3<ncell_sdr;j3++)
	{
	jj=xcell0[j3]+ycell0[j3]*nx*2;
	vtemp[jj+nstat+nevn] += lencell0[j3]*utemp[i]/sol[jj].wei;
#if ANIWEI
			vtemp[jj+nstat+nevn+nx*ny*4] += lencell0[j3]*utemp[i]*cos(2.0*PIO180*phase[i].azse)/sol[jj].weicos;
			vtemp[jj+nstat+nevn+2*nx*ny*4] += lencell0[j3]*utemp[i]*sin(2.0*PIO180*phase[i].azse)/sol[jj].weisin;
#endif
	}
	vtemp[phase[i].stnno] += utemp[i]/station[phase[i].stnno].wei;
	vtemp[phase[i].evnno+nstat] += utemp[i]/event[phase[i].evnno].wei;
	ii++;
}
xmesh=pow(2,grid_time)*xmesh0/2.0;
ymesh=pow(2,grid_time)*ymesh0/2.0;
nx=(xmax-xmin)*xmesh;
ny=(ymax-ymin)*ymesh;
//fprintf(stderr,"add======================grid over\n");
//========================================================================================================		
	/*zero station delays*/
	for(i=0;i<nstat;i++) {
		if(station[i].count==0) continue;
#if ANIWEI
		vtemp[i] += utemp[npha+3*nx*ny+1]/station[i].wei;
#else
		vtemp[i] += utemp[npha+nx*ny+1]/station[i].wei;
#endif
	}
	
	
	iii=ii;
	
	/*damp slowness here*/
	/*divide by weight*/
	for(ii=0,j=0;j<ny;j++) for(i=0;i<nx;i++) {
		if(sol[i+nx*j].nu==0) continue;
		vtemp[i+nx*j+nstat+nevn] += -4.0*weight*utemp[npha+i+nx*j]/sol[i+nx*j].wei;
		ii++;		
//--------------------------------------------------------------------------------------
		jj2=i/2+(nx/2)*(j/2);
		jj22=(i/2)*2+nx*(j/2)*2;	
		//=================================================1/2
		if(raysnum[0][jj2]<numax&&raysnum[0][jj2]>=0)
		{
		//left	
		if(i-2>=0 && sol[i+nx*j-1].nu!=0) 
		{
		vtemp[nstat+nevn+i+nx*j-1] += 0.5*weight*utemp[npha+i+nx*j]/sol[i+nx*j-1].wei;
		vtemp[nstat+nevn+i+nx*(j+1)-1] += 0.5*weight*utemp[npha+i+nx*j]/sol[i+nx*(j+1)-1].wei;
		} 
		if(i-2>=0 && sol[i+nx*j-1].nu==0&& sol[i+nx*j-2].nu!=0) vtemp[nstat+nevn+i+nx*j-2] += weight*utemp[npha+i+nx*j]/sol[i+nx*j-2].wei;
		//right
		if(i+2<nx && raysnum[0][jj2+1]>=numax) 
		{
		vtemp[nstat+nevn+i+nx*j+2] += 0.5*weight*utemp[npha+i+nx*j]/sol[i+nx*j+2].wei;
		vtemp[nstat+nevn+i+nx*(j+1)+2] += 0.5*weight*utemp[npha+i+nx*j]/sol[i+nx*(j+1)+2].wei;
		}
		if(i+2<nx && raysnum[0][jj2+1]>0&&raysnum[0][jj2+1]<numax) vtemp[nstat+nevn+i+nx*j+2] += weight*utemp[npha+i+nx*j]/sol[i+nx*j+2].wei;
		//down	
		if(j-2>=0 && sol[i+nx*(j-1)].nu!=0) 
		{
		vtemp[nstat+nevn+i+nx*(j-1)] += 0.5*weight*utemp[npha+i+nx*j]/sol[i+nx*(j-1)].wei;
		vtemp[nstat+nevn+i+nx*(j-1)+1] += 0.5*weight*utemp[npha+i+nx*j]/sol[i+nx*(j-1)+1].wei;
		}
		if(j-2>=0 && sol[i+nx*(j-1)].nu==0&& sol[i+nx*(j-2)].nu!=0) vtemp[nstat+nevn+i+nx*(j-2)] += weight*utemp[npha+i+nx*j]/sol[i+nx*(j-2)].wei;	
		//up	
		if(j+2<ny && raysnum[0][jj2+nx/2]>=numax) 
		{
		vtemp[nstat+nevn+i+nx*(j+2)] += 0.5*weight*utemp[npha+i+nx*j]/sol[i+nx*(j+2)].wei;
		vtemp[nstat+nevn+i+nx*(j+2)+1] += 0.5*weight*utemp[npha+i+nx*j]/sol[i+nx*(j+2)+1].wei;
		} 
		if(j+2<ny && raysnum[0][jj2+nx/2]>=0 && raysnum[0][jj2+nx/2]<numax) vtemp[nstat+nevn+i+nx*(j+2)] += weight*utemp[npha+i+nx*j]/sol[i+nx*(j+2)].wei;
		}
		
//=================================================1/4
		if(raysnum[0][jj2]>=numax)		
		{
		//left		
		if(i-1>=0 && sol[i+nx*j-1].nu!=0) vtemp[i+nx*j-1+nstat+nevn] += weight*utemp[npha+i+nx*j]/sol[i+nx*j-1].wei;
		if(i-2>=0 && sol[i+nx*j-1].nu==0 && sol[jj22-2].nu!=0) vtemp[jj22-2+nstat+nevn] += weight*utemp[npha+i+nx*j]/sol[jj22-2].wei;
		//right
		if(i+1<nx && sol[i+nx*j+1].nu!=0) vtemp[i+nx*j+1+nstat+nevn] += weight*utemp[npha+i+nx*j]/sol[i+nx*j+1].wei;
		if(i+2<nx && sol[i+nx*j+1].nu==0 && sol[jj22+2].nu!=0) vtemp[jj22+2+nstat+nevn] += weight*utemp[npha+i+nx*j]/sol[jj22+2].wei;
		//down	
		if(j-1>=0 )
		{
			if(sol[i+nx*(j-1)].nu!=0) vtemp[i+nx*(j-1)+nstat+nevn] += weight*utemp[npha+i+nx*j]/sol[i+nx*(j-1)].wei;
			if(j>=2) {if ( sol[i+nx*(j-1)].nu==0&& sol[jj22-2*nx].nu!=0) vtemp[jj22-2*nx+nstat+nevn] += weight*utemp[npha+i+nx*j]/sol[jj22-2*nx].wei;}
		}
		//up
		if(j+1<ny && sol[i+nx*(j+1)].nu!=0) vtemp[i+nx*(j+1)+nstat+nevn] += weight*utemp[npha+i+nx*j]/sol[i+nx*(j+1)].wei;
		if(j+2<ny && sol[i+nx*(j+1)].nu==0 &&  sol[jj22+2*nx].nu!=0) vtemp[jj22+2*nx+nstat+nevn] += weight*utemp[npha+i+nx*j]/sol[jj22+2*nx].wei;
		}
	}

#if ANIWEI
	/*damp cos here*/
	for(ii=0,j=0;j<ny;j++) for(i=0;i<nx;i++) {
		if(sol[i+nx*j].nu==0) continue;
		vtemp[i+nx*j+nstat+nevn+nx*ny] += -4.0*aniwei*utemp[npha+i+nx*j+nx*ny]/sol[i+nx*j].weicos;
		ii++;
//--------------------------------------------------------------------------------------		
		jj2=i/2+(nx/2)*(j/2);
		jj22=(i/2)*2+nx*(j/2)*2;
		//=================================================1/2
		if(raysnum[0][jj2]<numax&&raysnum[0][jj2]>=0)
		{
		//left	
		if(i-2>=0 && sol[i+nx*j-1].nu!=0) 
		{
		vtemp[nstat+nevn+nx*ny+i+nx*j-1] += 0.5*aniwei*utemp[npha+i+nx*j+nx*ny]/sol[i+nx*j-1].weicos;
		vtemp[nstat+nevn+nx*ny+i+nx*(j+1)-1] += 0.5*aniwei*utemp[npha+i+nx*j+nx*ny]/sol[i+nx*(j+1)-1].weicos;
		} 
		if(i-2>=0 && sol[i+nx*j-1].nu==0&& sol[i+nx*j-2].nu!=0) vtemp[nstat+nevn+nx*ny+i+nx*j-2] += aniwei*utemp[npha+i+nx*j+nx*ny]/sol[i+nx*j-2].weicos;
		//right
		if(i+2<nx && raysnum[0][jj2+1]>=numax) 
		{
		vtemp[nstat+nevn+nx*ny+i+nx*j+2] += 0.5*aniwei*utemp[npha+i+nx*j+nx*ny]/sol[i+nx*j+2].weicos;
		vtemp[nstat+nevn+nx*ny+i+nx*(j+1)+2] += 0.5*aniwei*utemp[npha+i+nx*j+nx*ny]/sol[i+nx*(j+1)+2].weicos;
		}
		if(i+2<nx && raysnum[0][jj2+1]>0&&raysnum[0][jj2+1]<numax) vtemp[nstat+nevn+nx*ny+i+nx*j+2] += aniwei*utemp[npha+i+nx*j+nx*ny]/sol[i+nx*j+2].weicos;
		//down	
		if(j-2>=0 && sol[i+nx*(j-1)].nu!=0) 
		{
		vtemp[nstat+nevn+nx*ny+i+nx*(j-1)] += 0.5*aniwei*utemp[npha+i+nx*j+nx*ny]/sol[i+nx*(j-1)].weicos;
		vtemp[nstat+nevn+nx*ny+i+nx*(j-1)+1] += 0.5*aniwei*utemp[npha+i+nx*j+nx*ny]/sol[i+nx*(j-1)+1].weicos;
		}
		if(j-2>=0 && sol[i+nx*(j-1)].nu==0&& sol[i+nx*(j-2)].nu!=0) vtemp[nstat+nevn+nx*ny+i+nx*(j-2)] += aniwei*utemp[npha+i+nx*j+nx*ny]/sol[i+nx*(j-2)].weicos;	
		//up	
		if(j+2<ny && raysnum[0][jj2+nx/2]>=numax) 
		{
		vtemp[nstat+nevn+nx*ny+i+nx*(j+2)] += 0.5*aniwei*utemp[npha+i+nx*j+nx*ny]/sol[i+nx*(j+2)].weicos;
		vtemp[nstat+nevn+nx*ny+i+nx*(j+2)+1] += 0.5*aniwei*utemp[npha+i+nx*j+nx*ny]/sol[i+nx*(j+2)+1].weicos;
		} 
		if(j+2<ny && raysnum[0][jj2+nx/2]>=0 && raysnum[0][jj2+nx/2]<numax) vtemp[nstat+nevn+nx*ny+i+nx*(j+2)] += aniwei*utemp[npha+i+nx*j+nx*ny]/sol[i+nx*(j+2)].weicos;
		
		
		}
//=================================================1/4
		if(raysnum[0][jj2]>=numax)		
		{
		//left		
		if(i-1>=0 && sol[i+nx*j-1].nu!=0) vtemp[i+nx*j-1+nstat+nevn+nx*ny] += aniwei*utemp[npha+i+nx*j+nx*ny]/sol[i+nx*j-1].weicos;
		if(i-2>=0 && sol[i+nx*j-1].nu==0 && sol[jj22-2].nu!=0) vtemp[jj22-2+nstat+nevn+nx*ny] += aniwei*utemp[npha+i+nx*j+nx*ny]/sol[jj22-2].weicos;
		//right
		if(i+1<nx && sol[i+nx*j+1].nu!=0) vtemp[i+nx*j+1+nstat+nevn+nx*ny] += aniwei*utemp[npha+i+nx*j+nx*ny]/sol[i+nx*j+1].weicos;
		if(i+2<nx && sol[i+nx*j+1].nu==0 && sol[jj22+2].nu!=0) vtemp[jj22+2+nstat+nevn+nx*ny] += aniwei*utemp[npha+i+nx*j+nx*ny]/sol[jj22+2].weicos;
		//down	
		if(j-1>=0 )
		{
			if(sol[i+nx*(j-1)].nu!=0) vtemp[i+nx*(j-1)+nstat+nevn+nx*ny] += aniwei*utemp[npha+i+nx*j+nx*ny]/sol[i+nx*(j-1)].weicos;
			if(j>=2) {if ( sol[i+nx*(j-1)].nu==0&& sol[jj22-2*nx].nu!=0) vtemp[jj22-2*nx+nstat+nevn+nx*ny] += aniwei*utemp[npha+i+nx*j+nx*ny]/sol[jj22-2*nx].weicos;}
		}
		//up
		if(j+1<ny && sol[i+nx*(j+1)].nu!=0) vtemp[i+nx*(j+1)+nstat+nevn+nx*ny] += aniwei*utemp[npha+i+nx*j+nx*ny]/sol[i+nx*(j+1)].weicos;
		if(j+2<ny && sol[i+nx*(j+1)].nu==0 &&  sol[jj22+2*nx].nu!=0) vtemp[jj22+2*nx+nstat+nevn+nx*ny] += aniwei*utemp[npha+i+nx*j+nx*ny]/sol[jj22+2*nx].weicos;
		}
	}
	/*damp sin here*/
	for(ii=0,j=0;j<ny;j++) for(i=0;i<nx;i++) {
		if(sol[i+nx*j].nu==0) continue;
		vtemp[i+nx*j+nstat+nevn+2*nx*ny] += -4.0*aniwei*utemp[npha+i+nx*j+2*nx*ny]/sol[i+nx*j].weisin;
		ii++;

//--------------------------------------------------------------------------------------
		
		jj2=i/2+(nx/2)*(j/2);
		jj22=(i/2)*2+nx*(j/2)*2;
		//=================================================1/2
		if(raysnum[0][jj2]<numax&&raysnum[0][jj2]>=0)
		{
		//left	
		if(i-2>=0 && sol[i+nx*j-1].nu!=0) 
		{
		vtemp[nstat+nevn+2*nx*ny+i+nx*j-1] += 0.5*aniwei*utemp[npha+i+nx*j+2*nx*ny]/sol[i+nx*j-1].weisin;
		vtemp[nstat+nevn+2*nx*ny+i+nx*(j+1)-1] += 0.5*aniwei*utemp[npha+i+nx*j+2*nx*ny]/sol[i+nx*(j+1)-1].weisin;
		} 
		if(i-2>=0 && sol[i+nx*j-1].nu==0&& sol[i+nx*j-2].nu!=0) vtemp[nstat+nevn+2*nx*ny+i+nx*j-2] += aniwei*utemp[npha+i+nx*j+2*nx*ny]/sol[i+nx*j-2].weisin;
		//right
		if(i+2<nx && raysnum[0][jj2+1]>=numax) 
		{
		vtemp[nstat+nevn+2*nx*ny+i+nx*j+2] += 0.5*aniwei*utemp[npha+i+nx*j+2*nx*ny]/sol[i+nx*j+2].weisin;
		vtemp[nstat+nevn+2*nx*ny+i+nx*(j+1)+2] += 0.5*aniwei*utemp[npha+i+nx*j+2*nx*ny]/sol[i+nx*(j+1)+2].weisin;
		}
		if(i+2<nx && raysnum[0][jj2+1]>0&&raysnum[0][jj2+1]<numax) vtemp[nstat+nevn+2*nx*ny+i+nx*j+2] += aniwei*utemp[npha+i+nx*j+2*nx*ny]/sol[i+nx*j+2].weisin;
		//down	
		if(j-2>=0 && sol[i+nx*(j-1)].nu!=0) 
		{
		vtemp[nstat+nevn+2*nx*ny+i+nx*(j-1)] += 0.5*aniwei*utemp[npha+i+nx*j+2*nx*ny]/sol[i+nx*(j-1)].weisin;
		vtemp[nstat+nevn+2*nx*ny+i+nx*(j-1)+1] += 0.5*aniwei*utemp[npha+i+nx*j+2*nx*ny]/sol[i+nx*(j-1)+1].weisin;
		}
		if(j-2>=0 && sol[i+nx*(j-1)].nu==0&& sol[i+nx*(j-2)].nu!=0) vtemp[nstat+nevn+2*nx*ny+i+nx*(j-2)] += aniwei*utemp[npha+i+nx*j+2*nx*ny]/sol[i+nx*(j-2)].weisin;	
		//up	
		if(j+2<ny && raysnum[0][jj2+nx/2]>=numax) 
		{
		vtemp[nstat+nevn+2*nx*ny+i+nx*(j+2)] += 0.5*aniwei*utemp[npha+i+nx*j+2*nx*ny]/sol[i+nx*(j+2)].weisin;
		vtemp[nstat+nevn+2*nx*ny+i+nx*(j+2)+1] += 0.5*aniwei*utemp[npha+i+nx*j+2*nx*ny]/sol[i+nx*(j+2)+1].weisin;
		} 
		if(j+2<ny && raysnum[0][jj2+nx/2]>=0 && raysnum[0][jj2+nx/2]<numax) vtemp[nstat+nevn+2*nx*ny+i+nx*(j+2)] += aniwei*utemp[npha+i+nx*j+2*nx*ny]/sol[i+nx*(j+2)].weisin;
		
		
		}
//=================================================1/4
		if(raysnum[0][jj2]>=numax)	
		{
		//left		
		if(i-1>=0 && sol[i+nx*j-1].nu!=0) vtemp[i+nx*j-1+nstat+nevn+2*nx*ny] += aniwei*utemp[npha+i+nx*j+2*nx*ny]/sol[i+nx*j-1].weisin;
		if(i-2>=0 && sol[i+nx*j-1].nu==0 && sol[jj22-2].nu!=0) vtemp[jj22-2+nstat+nevn+2*nx*ny] += aniwei*utemp[npha+i+nx*j+2*nx*ny]/sol[jj22-2].weisin;
		//right
		if(i+1<nx && sol[i+nx*j+1].nu!=0) vtemp[i+nx*j+1+nstat+nevn+2*nx*ny] += aniwei*utemp[npha+i+nx*j+2*nx*ny]/sol[i+nx*j+1].weisin;
		if(i+2<nx && sol[i+nx*j+1].nu==0 && sol[jj22+2].nu!=0) vtemp[jj22+2+nstat+nevn+2*nx*ny] += aniwei*utemp[npha+i+nx*j+2*nx*ny]/sol[jj22+2].weisin;
		//down	
		if(j-1>=0 )
		{
			if(sol[i+nx*(j-1)].nu!=0) vtemp[i+nx*(j-1)+nstat+nevn+2*nx*ny] += aniwei*utemp[npha+i+nx*j+2*nx*ny]/sol[i+nx*(j-1)].weisin;
			if(j>=2) {if ( sol[i+nx*(j-1)].nu==0&& sol[jj22-2*nx].nu!=0) vtemp[jj22-2*nx+nstat+nevn+2*nx*ny] += aniwei*utemp[npha+i+nx*j+2*nx*ny]/sol[jj22-2*nx].weisin;}
		}
		//up
		if(j+1<ny && sol[i+nx*(j+1)].nu!=0) vtemp[i+nx*(j+1)+nstat+nevn+2*nx*ny] += aniwei*utemp[npha+i+nx*j+2*nx*ny]/sol[i+nx*(j+1)].weisin;
		if(j+2<ny && sol[i+nx*(j+1)].nu==0 &&  sol[jj22+2*nx].nu!=0) vtemp[jj22+2*nx+nstat+nevn+2*nx*ny] += aniwei*utemp[npha+i+nx*j+2*nx*ny]/sol[jj22+2*nx].weisin;
		}
	}
#endif
	
	for(i=0;i<grid_time;i++)
	free(raysnum[i]);
	free(raysnum);
	/*fprintf(stderr,"atupv done\n");*/
}

/*************************************************/
double ran1(idum)
int idum;
/*returns random number 0.0 inclusive to 1.0 exclusive*/
/*thearn 6/90*/
{
	static int seed=0;
	float x;
	
	if(seed==0) {
		srand(idum);
		seed=1;
	}

	x = (float)rand()/(2147483647);
	//x = (float)rand()/(2147483647);	
	/*max on sun is 2**31-1*/
	/*max under gcc is 32767 */
	/*max under linux is 2147483647*/
	/*fprintf(stderr,"ran1:x,seed,idum %f %d %d\n",x,seed,idum);*/
	return x;
}

/***************************************************************/
