#include<stdio.h>
#include<string.h>

/*
c NMAX  = maximum number of bodies
c CMAX  = maximum number of close-encounter minima monitored simultaneously
c NMESS = maximum number of messages in message.in
c HUGE  = an implausibly large number
c NFILES = maximum number of files that can be open at the same time
*/
//      integer NMAX, CMAX, NMESS, NFILES
//      real*8 HUGE

#define NMAX 2000
#define CMAX 50
#define NMESS 200
#define HUGE 9.9e29
#define NFILES 50
/*
c DR = conversion factor from degrees to radians
c K2 = Gaussian gravitational constant squared
c AU = astronomical unit in cm
c MSUN = mass of the Sun in g
*/

#define PI 3.141592653589793e0
#define TWOPI (PI * 2.e0)
#define PIBY2 (PI * .5e0)
#define DR (PI / 180.e0)
#define K2 2.959122082855911e-4
#define AU 1.4959787e13
#define MSUN 1.9891e33

#define FUNSIG_MDT_HY (double* time, double* tstart, double* h0, double* tol, double* rmax, double en[], double am[], double jcen[], double* rcen, int* nbod, int* nbig, double m[], double x[][3], double v[][3], double s[][3], double rphys[], double rcrit[], double rce[3], int stat, char id[], double ngf, int *algor, int opt[], int* dtflag, int* ngflag, int* opflag, int* colflag, int* nclo, int iclo[], int jclo[], double dclo[], double tclo[], double ixvclo[][6], double jxvclo[][6], char outfile[][80], char mem[][80], int lmem[])

#define FUNSIG_MCO_H2DH (double* time, double jcen[], int* nbod, int* nbig, double* h, double m[], double xh[][3], double vh[][3], double x[][3], double v[][3], double ngf[][4], int* ngflag, int opt[])

#define FUNSIG_MCO_DH2H (double* time, double jcen[], int* nbod, int* nbig, double* h, double m[], double x[][3], double v[][3], double xh[][3], double vh[][3], double ngf[][4], int* ngflag, int opt[])

extern void mdt_hy_(double* time, double* tstart, double* h0, double* tol, double* rmax, double en[], double am[], double jcen[], double* rcen, int* nbod, int* nbig, double m[], double x[][3], double v[][3], double s[][3], double rphys[], double rcrit[], double rce[3], int stat, char id[], double ngf, int *algor, int opt[], int* dtflag, int* ngflag, int* opflag, int* colflag, int* nclo, int iclo[], int jclo[], double dclo[], double tclo[], double ixvclo[][6], double jxvclo[][6], char outfile[][80], char mem[][80], int lmem[]);

extern void mco_h2dh_(double* time, double jcen[], int* nbod, int* nbig, double* h, double m[], double xh[][3], double vh[][3], double x[][3], double v[][3], double ngf[][4], int* ngflag, int opt[]);

extern void mco_dh2h_(double* time, double jcen[], int* nbod, int* nbig, double* h, double m[], double x[][3], double v[][3], double xh[][3], double vh[][3], double ngf[][4], int* ngflag, int opt[]);

extern void mio_in_(double *time, double *tstart, double *tstop, double *dtout, int *algor, double *h0, double *tol, double *rmax, double *rcen, double jcen[], double en[], double am[], double *cefac, int *ndump, int *nfun, int *nbod, int *nbig, double m[], double x[][3], double v[][3], double s[][3], double rho[], double rceh[], int stat[], char id[][8], double epoch[], double ngf[][4], int opt[], int *opflag, int *ngflag, char outfile[][80], char dumpfile[][80], int lmem[], char mem[][80]);

extern void mxx_sync_(double* time, double* tstart, double* h0, double* tol, double jcen[], int* nbod, int* nbig, double m[], double x[][3], double v[][3], double s[][3], double rho[], double rceh[], int stat[], char id[][8], double epoch[], double ngf[][4], int opt[], int *ngflag);

extern void mal_hcon_(double *time, double *tstart, double *tstop, double *dtout, int *algor, double *h0, double *tol, double jcen[], double *rcen, double *rmax, double en[], double am[], double *cefac, int *ndump, int *nfun, int *nbod, int *nbig, double m[], double x[][3], double v[][3], double s[][3], double rho[], double rceh[], int stat[], char id[][8], double ngf[][4], int opt[], int *opflag, int *ngflag, char outfile[][80], char dumpfile[][80], char mem[][80], int lmem[], void (*mdt_hy_) FUNSIG_MDT_HY, void (*mco_h2dh_) FUNSIG_MCO_H2DH, void (*mco_dh2h) FUNSIG_MCO_DH2H);

extern void mio_dump_(double *time, double *tstart, double *tstop, double *dtout, int *algor, double *h0, double *tol, double jcen[], double *rcen, double *rmax, double en[], double am[], double *cefac, int *ndump, int *nfun, int *nbod, int *nbig, double m[], double x[][3], double v[][3], double s[][3], double rho[], double rceh[], int stat[], char id[][8], double ngf[][4], double epoch[], int opt[], int *opflag, char dumpfile[][80], char mem[][80], int lmem[]);

extern void mxx_en_(double jcen[], int* nbod, int* nbig, double m[], double xh[][3], double vh[][3], double s[][3], double* e, double* l2);

int main()
{
	int j,algor,nbod,nbig, stat[NMAX], lmem[NMESS];
	int opflag,ngflag,ndump,nfun;
	double m[NMAX],xh[NMAX][3],vh[NMAX][3],s[NMAX][3],rho[NMAX];
	double rceh[NMAX],epoch[NMAX],ngf[NMAX][4],rmax,rcen,jcen[3];
	double cefac,time,tstart,tstop,dtout,h0,tol,en[3],am[3];
	char id[NMAX][8];
	char outfile[3][80], dumpfile[4][80], mem[NMESS][80];
//      external mdt_mvs, mdt_bs1, mdt_bs2, mdt_ra15, mdt_hy, mdt_cb
//      external mdt_wb,mco_dh2h,mco_h2dh,mco_wb2h,mco_h2wb,mco_cb2h
//      external mco_h2cb,mco_b2h,mco_h2b,mco_h2mvs,mco_mvs2h,mco_iden

	int opt[8] = {0,1,1,2,0,1,0,0};

	mio_in_(&time, &tstart, &tstop, &dtout, &algor, &h0, &tol, &rmax, &rcen, jcen, en, am, &cefac, &ndump, &nfun, &nbod, &nbig, m, xh, vh, s, rho, rceh, stat, id, epoch, ngf, opt, &opflag, &ngflag, outfile, dumpfile, lmem, mem);


	if(opflag==-2)
	{
		//Ignoring some output
		mxx_sync_(&time, &tstart, &h0, &tol, jcen, &nbod, &nbig, m, xh, vh, s, rho, rceh, stat, id, epoch, ngf, opt, &ngflag);
		opflag = -1;
	}

	mal_hcon_(&time, &tstart, &tstop, &dtout, &algor, &h0, &tol, jcen, &rcen, &rmax, en, am, &cefac, &ndump, &nfun, &nbod, &nbig, m, xh, vh, s, rho, rceh, stat, id, ngf, opt, &opflag, &ngflag, outfile, dumpfile, mem, lmem, &mdt_hy_, &mco_h2dh_, &mco_dh2h_);

	for(j=2; j<=nbod; j++)
		epoch[j] = time;

	mio_dump_(&time, &tstart, &tstop, &dtout, &algor, &h0, &tol, jcen, &rcen, &rmax, en, am, &cefac, &ndump, &nfun, &nbod, &nbig, m, xh, vh, s, rho, rceh, stat, id, ngf, epoch, opt, &opflag, dumpfile, mem, lmem);

	//Ignoring some output
	mxx_en_(jcen, &nbod, &nbig, m, xh, vh, s, &(en[2]), &(am[2]));

	//Ignoring some output

	return 0;
}


	/*
		printf("--->%d\n", strlen(outfile[0]));

	 */
	//sprintf("----->%s\n", outfile[0]);
	/*
		int i;
		for(i=0; i<10; i++)
		{
		printf("%c", outfile[0][i]);
		if(outfile[0][i]=='\n') break;

		}
		printf("--->%d\n", i);
		printf("---->algor=%d opflag=%d\n", algor, opflag);
	 */
