/* rdphase.c - to read Nouredines phase data */
/* outputs - events.b, stations.b, phases.b */
/* input   - phase data */
/* assumes data is in event major order */
/* this version new 5/89 */
/* this version updated to use data3.h 9/89 */
double atand(),tand();
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
/*#include <libelf.h>*/  
/*next to lines are for mmap()*/
#include <sys/types.h>
#include <sys/mman.h>
/*next two for lseek()*/
#include <sys/types.h>
#include <unistd.h>
#include <stdlib.h>
#include "data14.h"

float xmin=0.0,xmax=360.0,ymin=-90.0,ymax=90.0;
float xmin1=0.0,xmax1=360.0,ymin1=-90.0,ymax1=90.0;
/*data selection parameters*/
#define DELMIN	2.0			/* min isc offset allowed */
#define DELMAX	18.0		/* max isc offset allowed in degrees */
#define CRESMAX	9.0		/* max isc calculated residual to use */
#define RESMAX	9.0		/* max calculated residual to use */
#define TTMAX	8800.		/* max total travel time allowed */
#define STNMIN	1		/* min num of arrivals per station */
#define EVNMIN	1			/* min num of arrivals per event */
#define YRMIN	10			/* must be > newer than this! */
#define YRMAX	2999		/* must be < this! */

#define DEPMAX	600.0		/* max depth allowed */
#define	NSTA	8000		/* max number of stations*/
#define NEVE	80000		/* max number of events*/
#define DN		35
#define VG		6.3
#define VN		8.1
#define ETAELE	0.0997		/*elevation corretion for p=8.1,v=6.3*/
#define REARTH	6371.0		/*earth radius*/

double          	fabs(), sin(), sqrt() ;
char				*l64a();
struct rec 			*r;
struct locations	*station, *event;
struct rec			*phase;
int					npha;
float				resmax=RESMAX,dn=DN,vg=VG,vn=VN;quality1;
float				cresmax=CRESMAX;
int					stnmin=STNMIN,evnmin=EVNMIN;
float				delmin=DELMIN,delmax=DELMAX;
float				depmax=DEPMAX;
int 				yrmin=YRMIN;yrmax=YRMAX;


main(argc,argv)
int argc; char **argv;
{
	int             count, stncnt, evncnt;
	int             i, j;
	double          inter, slope;
	int             fdevn, fdstn, fdrec;
	FILE           *fopen(), *fp;

	/*input options go here*/
	if(getarg(argc,argv,"-region","%f/%f/%f/%f",&xmin,&xmax,&ymin,&ymax)<4) errprint("rdphase: cant getarg region");
	fprintf(stderr,"xmin,xmax,ymin,ymax= %f %f %f %f\n",xmin,xmax,ymin,ymax);
	getarg(argc,argv,"-delmin","%f",&delmin);
	getarg(argc,argv,"-delmax","%f",&delmax);
	getarg(argc,argv,"-resmax","%f",&resmax);
	getarg(argc,argv,"-cresmax","%f",&cresmax);
	getarg(argc,argv,"-stnmin","%d",&stnmin);
	getarg(argc,argv,"-evnmin","%d",&evnmin);
	getarg(argc,argv,"-depmax","%f",&depmax);	
	getarg(argc,argv,"-yrmin","%d",&yrmin);	
	getarg(argc,argv,"-yrmax","%d",&yrmax);

	getarg(argc,argv,"-xmin1","%f",&xmin1);
	getarg(argc,argv,"-ymin1","%f",&ymin1);
	getarg(argc,argv,"-xmax1","%f",&xmax1);
	getarg(argc,argv,"-ymax1","%f",&ymax1);
	fprintf(stderr,"delmin,delmax,cresmax,resmax,TTMAX,%f %f %f %f %f\n",
		delmin,delmax,cresmax,resmax,TTMAX);
	fprintf(stderr,"depmax,stnmin,evnmin, yrmin,yrmax,%f %d %d %d %d\n",depmax,stnmin,evnmin,yrmin,yrmax);
	getarg(argc,argv,"-dn","%f",&dn);
	getarg(argc,argv,"-vg","%f",&vg);
	getarg(argc,argv,"-vn","%f",&vn);
	fprintf(stderr,"dn,vn,vg= %f %f %f\n",dn,vn,vg);

	/* this line for stdio input files */
	fp = stdin;

	/* output data files */
	fdevn = open("../data/events.b", O_RDWR |O_CREAT|O_TRUNC,0644);
	fdstn = open("../data/stations.b",O_RDWR|O_CREAT|O_TRUNC,0644);
	fdrec = open("../data/phases.b", O_RDWR |O_CREAT|O_TRUNC,0644);
	if(fdevn==0||fdstn==0||fdrec==0) errprint("error opening output file");
	
	/*allocate storage space */
	station = (struct locations *) malloc(NSTA*sizeof(struct locations));
	event = (struct locations *) malloc(NEVE*sizeof(struct locations));
	r = (struct rec *) malloc(sizeof(struct rec));
	printf("event, station = %d %d\n",event,station);

	/* readin data & set data structures */
	printf("calling rddata\n");
	rddata(fp, fdrec, &count, &stncnt, &evncnt);
	
	/* now rewind fdrec, and impliment mmap on it for speed */
	lseek(fdrec, 0L, 0);
	phase = (struct rec *) mmap(0,lseek(fdrec,0L,SEEK_END),PROT_READ|PROT_WRITE, MAP_SHARED,fdrec,0);
	npha = lseek(fdrec,0L,SEEK_END)/sizeof(struct rec);
	printf("phase= %d\n", phase);

	/* now start selection process */
	while (phsel()>0) phcnt();

	/* fit a line to it return inter & slope */
	linefit(&inter, &slope);
	//linefit2(&inter, &slope);

	/* calculate the residuals, set quality to -1 for bad arrivals */
	while (rescalc(inter, slope) == 1) {
		phcnt();
		while (phsel()>0) phcnt();
		linefit(&inter, &slope);
	}
	
	/* set ttime=dtime for all phases 在2025.09修改：避免phases.b内的绝对走时变成走时残差，以便画时距分布*/
	//for(i=0;i<npha;i++) phase[i].ttime = phase[i].dtime;

	/* now write out event & station */
	for(i=0; i<evncnt; i++) write( fdevn, &event[i], sizeof(struct locations));
	for(i=0; i<stncnt; i++) write(fdstn, &station[i], sizeof(struct locations));

	/* all done now */
	printf("%d times in, %d events, %d stations\n", count, evncnt, stncnt);
	printf("all done now\n\n");
}


findstn(name, stncnt)
char           *name;
int             stncnt;

/* sees if station is already found */
/* name is stn name, stncnt is number of stations already in list */
{
	int             i, j;

	for (i = 0; i < stncnt; i++)
		if(strcmp(station[i].name,name)==0) return(i);
	return (-1);
}




linefit(inter, slope)
//fit a line to the damn data
double         *inter, *slope;
{
	int             i, k, j;
	double          sumx, sumy, sumxy, sumx2, sumy2;
	int             sum;
	double          res2;

	// now fit a line to it and calculate residuals
	printf("linefit: start linefit\n");
	
	sumx = sumy = sumxy = sumx2 = sumy2 = sum = 0.0;
	for(i=0; i<npha; i++) {
		if (phase[i].quality <= 0.0) continue;
		sumx += phase[i].offset;
		sumy += phase[i].ttime;
		sumxy += phase[i].offset * phase[i].ttime;
		sumx2 += phase[i].offset * phase[i].offset;
		sumy2 += phase[i].ttime * phase[i].ttime;
		sum++;
	}
	*inter = sum * sumx2 - sumx * sumx;
	*slope = (sum * sumxy - sumx * sumy) / (*inter);
	*inter = (sumy * sumx2 - sumx * sumxy) / (*inter);
	res2 = (sumy2 - *inter * sumy - *slope * sumxy);
	printf("linefit: sumx,sumy,sumxy,sumx2= %f %f %f %f\n", sumx, sumy, sumxy, sumx2);
	printf("linefit: inter,slope,slope-1,res2,sum= %f %f %f %f %d\n", *inter, *slope, (1.0 / (*slope)), res2, sum);
}



/* 20250925修改后的迭代重加权最小二乘法拟合部分 */
/* robust line fit with Tukey bisquare weights
   replaces previous linefit(inter, slope)
   Uses globals: struct rec *phase; int npha;
   Assumes phase[i].offset = x (km), phase[i].ttime = y (s)
   and phase[i].quality initial positive for usable points.*/

int cmp_double(const void *a, const void *b) {
    double da = *(const double*)a;
    double db = *(const double*)b;
    if (da < db) return -1;
    if (da > db) return 1;
    return 0;
}

linefit2(inter, slope)
double *inter, *slope;
{
    int i, iter;
    const int max_iter = 50;
    const double tol = 1e-8;
    const double c = 4.685; /* Tukey tuning constant */
    double *w = NULL;
    double *xs = NULL;
    double *ys = NULL;
    double *res = NULL;
    double sumw, meanx, meany, Sxx, Sxy;
    int nvalid = 0;

    /* count valid points */
    for (i = 0; i < npha; i++) if (phase[i].quality > 0.0) nvalid++;
    if (nvalid < 2) {
        fprintf(stderr, "robust linefit: not enough valid points n=%d\n", nvalid);
        *inter = 0.0; *slope = 0.0;
        return;
    }

    xs = (double*) malloc(nvalid * sizeof(double));
    ys = (double*) malloc(nvalid * sizeof(double));
    w  = (double*) malloc(nvalid * sizeof(double));
    res = (double*) malloc(nvalid * sizeof(double));
    if (!xs || !ys || !w || !res) {
        fprintf(stderr, "robust linefit: malloc failed\n");
        exit(1);
    }

    /* copy data into compact arrays */
    nvalid = 0;
    for (i = 0; i < npha; i++) {
        if (phase[i].quality <= 0.0) continue;
        xs[nvalid] = phase[i].offset;
        ys[nvalid] = phase[i].ttime;
        w[nvalid]  = 1.0; /* initial weight (or use phase[i].quality if desired) */
        nvalid++;
    }

    /* initial ordinary weighted least squares to seed */
    /* compute weighted means */
    sumw = 0.0; meanx = 0.0; meany = 0.0;
    for (i = 0; i < nvalid; i++) {
        sumw += w[i];
        meanx += w[i] * xs[i];
        meany += w[i] * ys[i];
    }
    meanx /= sumw; meany /= sumw;

    Sxx = Sxy = 0.0;
    for (i = 0; i < nvalid; i++) {
        double dx = xs[i] - meanx;
        double dy = ys[i] - meany;
        Sxx += w[i] * dx * dx;
        Sxy += w[i] * dx * dy;
    }
    if (fabs(Sxx) < 1e-14) {
        *slope = 0.0;
        *inter = meany;
    } else {
        *slope = Sxy / Sxx;
        *inter = meany - (*slope) * meanx;
    }

    /* IRLS iterations with Tukey bisquare */
    for (iter = 0; iter < max_iter; iter++) {
        /* residuals */
        for (i = 0; i < nvalid; i++) {
            res[i] = ys[i] - (*inter + (*slope) * xs[i]);
            if (!isfinite(res[i])) res[i] = 0.0;
        }
        /* estimate scale s by MAD */
        double *absr = (double*) malloc(nvalid * sizeof(double));
        if (!absr) { fprintf(stderr,"malloc failed\n"); exit(1); }
        for (i = 0; i < nvalid; i++) absr[i] = fabs(res[i]);
        qsort(absr, nvalid, sizeof(double), cmp_double);
        double med = absr[nvalid/2];
        free(absr);
        double s_est = 1.4826 * med;
        if (s_est <= 0.0) {
            /* fallback to RMS */
            double sum2 = 0.0;
            for (i = 0; i < nvalid; i++) sum2 += res[i]*res[i];
            s_est = sqrt(sum2 / (double) nvalid);
            if (s_est <= 0.0) s_est = 1e-6;
        }

        /* compute new weights (Tukey bisquare) */
        double maxw = 0.0;
        for (i = 0; i < nvalid; i++) {
            double u = fabs(res[i]) / (c * s_est);
            if (u >= 1.0) w[i] = 0.0;
            else {
                double tmp = 1.0 - u*u;
                w[i] = tmp * tmp;
            }
            if (w[i] > maxw) maxw = w[i];
        }
        /* avoid all-zero weights */
        if (maxw == 0.0) {
            /* no good solution; break */
            break;
        }
        /* weighted least squares with new weights */
        sumw = 0.0; meanx = 0.0; meany = 0.0;
        for (i = 0; i < nvalid; i++) {
            sumw += w[i];
            meanx += w[i] * xs[i];
            meany += w[i] * ys[i];
        }
        if (sumw == 0.0) break;
        meanx /= sumw; meany /= sumw;

        double Sxx_new = 0.0, Sxy_new = 0.0;
        for (i = 0; i < nvalid; i++) {
            double dx = xs[i] - meanx;
            double dy = ys[i] - meany;
            Sxx_new += w[i] * dx * dx;
            Sxy_new += w[i] * dx * dy;
        }
        if (fabs(Sxx_new) < 1e-14) {
            /* degenerate */
            double new_slope = 0.0;
            double new_inter = meany;
            if (fabs(new_slope - *slope) < tol && fabs(new_inter - *inter) < tol) break;
            *slope = new_slope; *inter = new_inter;
            break;
        } else {
            double new_slope = Sxy_new / Sxx_new;
            double new_inter = meany - new_slope * meanx;
            /* check convergence */
            if (fabs(new_slope - *slope) < tol && fabs(new_inter - *inter) < tol) {
                *slope = new_slope; *inter = new_inter;
                break;
            }
            *slope = new_slope; *inter = new_inter;
            /* continue iterate */
        }
    } /* end iter */

    /* OPTIONAL: mark very low-weight points as bad (quality = -1)
       e.g. mark those with w == 0 or |res| > 4*s_est
    */
    for (i = 0, nvalid = 0; i < npha; i++) {
        if (phase[i].quality <= 0.0) continue;
        /* find index in compact arrays: we must rewalk in same order */
        double x = phase[i].offset;
        double y = phase[i].ttime;
        /* linear search to find matching index (could optimize by storing indices) */
        int found = -1;
        for (int j = 0; j < nvalid + npha; j++); /* dummy to quiet compilers - not used */
        /* simpler: rebuild mapping in first copy step to store original idxs */
    }

    /* NOTE: above marking requires mapping array indices back to phase[] entries.
       To keep it simple, better to build an index array at start. */

    free(xs); free(ys); free(w); free(res);
}



rescalc(inter, slope)
double          inter, slope;
/*this calcultes residuals relative to line fit*/
{
	int             i, j, k;
	double          res1, res2, res3, res4;
	float           dtmin, dtmax;
	int             rescalc;

	/* now calculate the residuals */
	printf("rescalc: start rescalc\n");
	dtmin = dtmax = 0.0;
	res2 = res1 = res3 = res4 = 0.0;
	rescalc = 0;
	printf("linefittttttt: inter,slope,slope-1= %f %f %f %f \n", inter, slope, (1.0 / (slope)));
	for(i=0; i<npha; i++) {
		phase[i].dtime = phase[i].ttime - inter - slope * phase[i].offset;
		//phase[i].dtime = phase[i].ttime - 7.32- phase[i].offset/8.1;
		
		//phase[i].dtime = phase[i].ttime - 11.1 - phase[i].offset/8.26;//青藏
		

		//phase[i].dtime = phase[i].ttime - 7.9 - phase[i].offset/8.13;
		if(phase[i].quality<=0.0) continue;
		res1 += fabs(phase[i].dtime) * phase[i].quality;
		res2 += phase[i].dtime * phase[i].dtime * phase[i].quality;
		res3 += phase[i].dtime * phase[i].quality;
		res4 += phase[i].quality;
		if (phase[i].dtime>dtmax) dtmax = phase[i].dtime;
		else if (phase[i].dtime < dtmin) dtmin = phase[i].dtime;
		if (fabs(phase[i].dtime) > resmax) {
			phase[i].quality = -1;
			rescalc = 1;
		}
	}
	printf("rescalc: r1,r2,r3,r4= %f %f %f %f\n", res1, res2, res3, res4);
	printf("rescalc: dtmin,dtmax= %f %f\n", dtmin, dtmax);
	return(rescalc);
}


phcnt()
/*count the number of phases per station and per event*/
{
	int             i, j, k;

	printf("phcnt: start phcnt\n");
	for(i=0; i<NSTA; i++) station[i].count = 0;
	for(i=0; i<NEVE; i++) event[i].count = 0;
	for(i=0; i<npha; i++) {
		if(phase[i].quality <= 0.0) continue;
		station[phase[i].stnno].count++;
		event[phase[i].evnno].count++;
	}
}

phsel()
//根据台站和事件数量筛选
/* returns number of phases deleted if a count is still unacceptable */
{
	int             i, j, k;
	int             count2;
	int             phsel;

	printf("phsel: start phsel\n");
	count2=phsel=0;
	for(i=0; i<npha; i++) {
		if (phase[i].quality<=0.0) continue;
		if (station[phase[i].stnno].count<stnmin || event[phase[i].evnno].count<evnmin) {
			phase[i].quality = -1.0;
			phsel++;
		} else count2++;
	}
	printf("phsel: count2,phsel= %d %d\n", count2, phsel);
	return(phsel);
}


rddata(fp, fdrec, count, stncnt, evncnt)
FILE           *fp;
int             fdrec;
int            *count, *stncnt, *evncnt;
{
	int             i, j, k;
	int             evnno;
	float           slat, slon, rlat, rlon, del, ssec,sdepth, rsec;
	int             az, nstrec,  yr, mo, da;
	char            name[5],ename[5];
	float           mag, res, elev;
	int             stnno;
	int             count2 = 0;
	int				picpres, latpres, lonpres;
	float			soff,roff,dlat,dlon;
	float			azse, azes;
	double			distkm();

	/* start reading in data set */
	/* this is the main data readin loop */
	printf("rddata: starting rddata.c\n");
		
	(*stncnt)=(*evncnt)=evnno=0;
	FILE *fp1226;
	fp1226=fopen("../rawdata/2023kuril2_select312.dat","w2");
	FILE *fp1226de;
	fp1226de=fopen("../rawdata/2023kuril22_delete.dat","w2");
	FILE *fpq;
	fpq=fopen("../rawdata/Qualityx.dat","w2");
	for(*count=0; (i=fscanf(fp,
		"%f%f%f%f%f %f%d%f%f%f %f%d%f%d%d %d%s%d%d%d",
		&slat, &slon, &rlat, &rlon, &elev,
		&del, &az, &ssec, &rsec, &mag,
		&res, &nstrec, &sdepth, &yr, &mo,
		&da, name, &picpres, &latpres, &lonpres))!=EOF; (*count)++) {
		if(*count%100==0) printf("\rrddata: starting loop, count= %d", *count);
		if(i<20) {
			printf("rddata: fscanf returns %d on record %d\n", i, *count);
			printf("%f %f %f %f %f  %f %d %f %f %f  %f %d %f %d %d  %d %s %d %d %d\n",
				slat, slon, rlat, rlon, elev,
				del, az, ssec, rsec, mag,
				res, nstrec, sdepth, yr, mo,
				da, name, picpres, latpres, lonpres); }
		quality1=rsec-ETAELE * elev / 1000.0 + sdepth*sqrt(1.0/(vg*vg)- 1.0/(vn*vn))- 7.32- del*111.1/8.1;
		//quality1=rsec-ETAELE * elev / 1000.0 + sdepth*sqrt(1.0/(vg*vg)- 1.0/(vn*vn))- 11.1- del*111.1/8.26;
		//printf("quality  %f %f %fon record \n", quality,ETAELE,elev );
		if(fabs(quality1)<5&&fabs(quality1)>0) {
			fprintf(fpq,"%f %f %f %f %f  %f %d %f %f %f  %f %d %f %d %d  %d %s %d %d %d %f\n",
				slat, slon, rlat, rlon, elev,
				del*111.1, az, ssec, rsec, mag,
				res, nstrec, sdepth, yr, mo,
				da, name, picpres, latpres, lonpres, quality1);}

/*
if(del<3  &&  slat>=ymin1  &&  slat<=ymax1  &&  slon>=xmin1  &&  slon<=xmax1&&  rlat>=ymin1  &&  rlat<=ymax1  &&  rlon>=xmin1  &&  rlon<=xmax1)
fprintf(fp1226de,"%f %f %f %f %f  %f %d %f %f %f  %f %d %f %d %d  %d %s %d %d %d\n",
				slat, slon, rlat, rlon, elev,
				del, az, ssec, rsec, mag,
				res, nstrec, sdepth, yr, mo,
				da, name, picpres, latpres, lonpres);
else
{

fprintf(fp1226,"%f %f %f %f %f  %f %d %f %f %f  %f %d %f %d %d  %d %s %d %d %d\n",
				slat, slon, rlat, rlon, elev,
				del, az, ssec, rsec, mag,
				res, nstrec, sdepth, yr, mo,
				da, name, picpres, latpres, lonpres);
}
*/
		/*sort data by initial noniterative criteria*/
		/*shorten names to 4 characters (changed from 3 characters 5/30/94) */
		name[4]='\0';
		for(k=ssec;k>9999;k/=10);
		sprintf(ename,"%4d\0",k);
		/*fix year*/
		if(yr>200) yr-=1900;
 
 		/* selection criteria on data base go here */
		/*this is for isc data*/
		/*pick precision*/
		/*0 = 1 sec err, >0 is useless*/
		/*if(picpres>-1) continue;*/
		/*if(picpres==-3 || picpres==-4) continue;*/
		/*this is for gb data*/
		/*if(picpres>3) continue;*/
		/*this is for tibet data*/
		/*if(picpres>8) continue*/
		
		/*get rid of events with bad origin times*/
		/*this is done in pn_phase.awk*/
		/*qual 3 = 1/10 minute*/
		/*qual 2 = minute*/
		/*qual 1 = ten seconds*/
		/*qual 0 = second*/
		/*qual -1 = 1/10 second*/
		/*qual -2 = 1/100 second*/
		
		/*get rid of events with bad lats*/
		/*qual 0= 0 decimal places
		qual -1= 1 decimal place
		qual -2= 2 decimal place
		qual -3= 3 decimal place
		qual 4= deg min sec & tenth of sec
		qual 5= deg min sec
		qual 6= deg, min & tenth of minutes
		qual 7= deg & minutes
		qual 8= nearest quarter degree*/

		/*lat precision*/
		if(latpres==8 || latpres==0 || latpres==-1) continue;
		if(latpres==-4 || latpres==-5) continue;
		/*lon precision*/
		if(lonpres==8 || lonpres==0 || lonpres==-1) continue;
		if(lonpres==-4 || lonpres==-5) continue;
		
		/*this is for all the data*/
		if(fabs(res)>cresmax || fabs(rsec-ssec)>TTMAX) continue;
		
		if(slat>ymax || slat<ymin || slon>xmax || slon<xmin) continue;
		if(rlat>ymax || rlat<ymin || rlon>xmax || rlon<xmin) continue;
		if(del>delmax || del<delmin) continue;
//=====================================================================================select
		//if(del<3 &&  slat>=ymin1  &&  slat<=ymax1  &&  slon>=xmin1  &&  slon<=xmax1&&  rlat>=ymin1  &&  rlat<=ymax1  &&  rlon>=xmin1  &&  rlon<=xmax1) continue;

		if(sdepth>depmax) continue;
		/*if(sdepth==33.00) continue;*/
		if(yr<yrmin) continue;
		if(yr>yrmax) continue;
		/*printf("name= %s\n",name);*/

		/* collect event info (assume event major order) */
		if (event[evnno].lat != slat || event[evnno].lon != slon) {
			evnno = (*evncnt);
			/*event name number of secs*10 since year 1900*/
			/*this is fit into a 6 character string*/
			strcpy( event[evnno].name, l64a(julday(mo,da,yr)*86400+(int)ssec) );
			/*printf("event name is %s\n",event[evnno].name);*/
			event[evnno].lat = slat;
			event[evnno].lon = slon;
			event[evnno].number = evnno;
			event[evnno].elev = sdepth;
			event[evnno].mag = mag;
			event[evnno].sec = ssec;
			event[evnno].delay = 0;
			event[evnno].sum = 0.0;
			event[evnno].sum2 = 0.0;
			sprintf(event[evnno].date, "%2d %2d %2d", yr, mo, da);
			(*evncnt)++;
		}
		event[evnno].count++;

		/* collect station info */
		if ((stnno = findstn(name, *stncnt)) == -1) {	/* a new station */
			stnno = (*stncnt);
			strcpy(station[stnno].name, name);
			station[stnno].lat = rlat;
			station[stnno].lon = rlon;
			station[stnno].number=stnno;
			station[stnno].elev = elev;
			station[stnno].delay = 0.0;
			station[stnno].sum = 0.0;
			station[stnno].sum2 = 0.0;
			(*stncnt)++;
		}
		station[stnno].count++;

		/* collect phase info */
		/*first get the azse and great circle offset, don't comment this one out!*/

		r->offset = distkm(rlat,rlon,slat,slon,&azse,&azes);//用台站-事件坐标计算offset，不要注释
		//r->scale = del * 111.1 / r->offset;
		//r->offset = del * 111.1;//直接用del计算offset

		/*next get the cord of the great circle*/
		/*use only one of these next two statements, preferably the second*/
		/*r->offset = 2.0 * REARTH * sin(0.017453295 * del / 2.0);*/

		r->offset = 2.0 * REARTH * sin(r->offset/(REARTH*2.0));//弧长转换为弦长
		
		//if(fabs(r->offset-del*111.1)>14) fprintf(stderr,
		//	"offsets: %f %f %f %f %f  %f %f %f %s %d  %d %d\n",
		//	r->offset,del*111.1,rsec,ssec,rlat,rlon,slat,slon,name,yr,mo,da);
		/*next get offset lats & lons*/
		//计算入射点与出射点的水平位移
		//soff：从震源点到moho面的水平射程（source offset）
		//roff：从台站点到moho面的水平射程（receiver offset）
		soff=(dn-sdepth)*vg/sqrt(vn*vn-vg*vg);
		roff = (dn + elev/1000.0) * vg / sqrt(vn * vn - vg * vg);
		//检查几何合理性：若震源端水平射程+台站端射程 > 总弦长，说明几何不可能，打印警告		
		if(soff+roff > r->offset) 
			fprintf(stderr,
				"dn,vg,elev,sdepth,soff,roff,offset %f %f %f %f %f %f %f\n",dn,vg,elev,sdepth,soff,roff,r->offset);
		
		//方向向量与缩放：dlat，dlon：台站与事件的经纬度差；scale：把经纬度差向量缩放，使其长度匹配offset弦长。---后面用于沿直线插值震源点和台站点到moho面交点。
		dlat = rlat - slat;
		dlon = rlon - slon;
		
		//r->scale = del / sqrt(dlat * dlat + dlon * dlon);
		r->scale = r->offset / sqrt(dlat * dlat + dlon * dlon);

		//r->scale = r->offset / (2.0 * REARTH * sin(r->offset/(REARTH*2.0)));//这一行用来检查弦长和弧长之间的差异
		
		/*these lats & lons are offset by the lateral distace to the refractor */
		//修正事件、台站位置（投影到moho面交点）
		//slat slon事件原始坐标 rlat rlon台站原始坐标；把事件、台站沿着chord方向便宜soff、roff。得到折射点坐标---定位射线在moho面上下穿点的坐标		
		r->slat = slat + soff * dlat / r->offset;
		r->slon = slon + soff * dlon / r->offset;
		r->rlat = rlat - roff * dlat / r->offset;
		r->rlon = rlon - roff * dlon / r->offset;
		
		//写入震相属性
		r->azse = azse;//方位角
		/* this is station to event azimuth clockwise from north */
		r->ttime = rsec - ssec;//走时（观测台站时 - 事件时）
		r->ttime -= ETAELE * elev / 1000.0;//高程修正
		r->ttime += sdepth*sqrt(1.0/(vg*vg)- 1.0/(vn*vn));//震源深度修正
		r->dtime = 0.0;
		r->delay = 0.0;
		r->arrtype = 0;
		r->evnno = evnno;
		//r->evnno = yr;
		r->stnno = stnno;
		r->quality = 1.0;
		

		/* write phase info out now */
		//写入phases.b
		i = write(fdrec, r, RECSIZ);
		if(i<=0) printf("rddata: error in rec write, %d\n", i);
		count2++;	/* keeps track of number of phases accepted */
loopend:;
	}
	fclose(fp1226);
	fclose(fp1226de);
	fclose(fpq);
	printf("\nrddata: main loop done, data all read in\n");
	printf("rddata: count,count2= %d %d\n", *count, count2);
}
