#include<stdio.h>
#include<string.h>

//Making all parameters in mercury.inc as macros

/*
c Parameters that you may want to alter at some point:
c
c NMAX  = maximum number of bodies
c CMAX  = maximum number of close-encounter minima monitored simultaneously
c NMESS = maximum number of messages in message.in
c HUGE  = an implausibly large number
c NFILES = maximum number of files that can be open at the same time
*/

#define NMAX 2000
#define CMAX 50
#define NMESS 200
#define HUGE 9.9e29
#define NFILES 50

/*
c Constants:
c
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

//Macros for function signatures since function pointers are used as argument
#define FUNCSIG_MAL_HCON_1 (double* time, double* tstart, double* h0, double* tol, double* rmax, double en[], double am[], double jcen[], double* rcen, int* nbod, int* nbig, double m[], double x[][3], double v[][3], double s[][3], double rphys[], double rcrit[], double rce[3], int stat, char id[], double ngf, int *algor, int opt[], int* dtflag, int* ngflag, int* opflag, int* colflag, int* nclo, int iclo[], int jclo[], double dclo[], double tclo[], double ixvclo[][6], double jxvclo[][6], char outfile[][80], char mem[][80], int lmem[])

#define FUNCSIG_MAL_HCON_2 (double* time, double jcen[], int* nbod, int* nbig, double* h, double m[], double xh[][3], double vh[][3], double x[][3], double v[][3], double ngf[][4], int* ngflag, int opt[])

#define FUNCSIG_MAL_HCON_3 (double* time, double jcen[], int* nbod, int* nbig, double* h, double m[], double x[][3], double v[][3], double xh[][3], double vh[][3], double ngf[][4], int* ngflag, int opt[])

#define FUNCSIG_MAL_HVAR (double* time, double* h0, double* hdid, double* tol, double jcen[],int* nbod, int* nbig, double mass[], double x0[][3], double v0[][3], double s[][3], double rphys[], double rcrit[], double ngf[][4], int stat[], int dtflag, int ngflag, int opt[], int nce, int ice[], int jce[], void (*force) FUNCSIG_FORCE)

#define FUNCSIG_FORCE (double* time, double jcen[], int* nbod, int* nbig, double m[], double x[][3], double v[][3], double s[][3], double rcrit[], double a[][3], int stat[], double ngf[][4], int* ngflag, int opt[], int nce, int ice[], int jce[])


//Function definitions
//Functions with signature FUNCSIG_MAL_HCON_1
extern void mdt_mvs_ FUNCSIG_MAL_HCON_1;
extern void mdt_hy_ FUNCSIG_MAL_HCON_1;
extern void mdt_cb_ FUNCSIG_MAL_HCON_1;
extern void mdt_wb_ FUNCSIG_MAL_HCON_1;

//Functions with signature FUNCSIG_MAL_HCON_2
extern void mco_h2mvs_ FUNCSIG_MAL_HCON_2;
extern void mco_iden_ FUNCSIG_MAL_HCON_2;
extern void mco_h2dh_ FUNCSIG_MAL_HCON_2;
extern void mco_h2cb_ FUNCSIG_MAL_HCON_2;
extern void mco_h2wb_ FUNCSIG_MAL_HCON_2;

//Functions with signature FUNCSIG_MAL_HCON_3
extern void mco_mvs2h_ FUNCSIG_MAL_HCON_3;
extern void mco_dh2h_ FUNCSIG_MAL_HCON_3;
extern void mco_cb2h_ FUNCSIG_MAL_HCON_3;
extern void mco_wb2h_ FUNCSIG_MAL_HCON_3;

//Functions with signature FUNCSIG_MAL_HVAR
extern void mdt_bs1_ FUNCSIG_MAL_HVAR;
extern void mdt_bs2_ FUNCSIG_MAL_HVAR;
extern void mdt_ra15_ FUNCSIG_MAL_HVAR;


//Function which does input
extern void mio_in_(double *time, double *tstart, double *tstop, double *dtout, int *algor, double *h0, double *tol, double *rmax, double *rcen, double jcen[], double en[], double am[], double *cefac, int *ndump, int *nfun, int *nbod, int *nbig, double m[], double x[][3], double v[][3], double s[][3], double rho[], double rceh[], int stat[], char id[][8], double epoch[], double ngf[][4], int opt[], int *opflag, int *ngflag, char outfile[][80], char dumpfile[][80], int lmem[], char mem[][80]);


//Functions which do dynamics
extern void mal_hcon_(double *time, double *tstart, double *tstop, double *dtout, int *algor, double *h0, double *tol, double jcen[], double *rcen, double *rmax, double en[], double am[], double *cefac, int *ndump, int *nfun, int *nbod, int *nbig, double m[], double x[][3], double v[][3], double s[][3], double rho[], double rceh[], int stat[], char id[][8], double ngf[][4], int opt[], int *opflag, int *ngflag, char outfile[][80], char dumpfile[][80], char mem[][80], int lmem[], void (*mdt_hy_) FUNCSIG_MAL_HCON_1, void (*mco_h2dh_) FUNCSIG_MAL_HCON_2, void (*mco_dh2h) FUNCSIG_MAL_HCON_3);

extern void mal_hvar_(double *time, double *tstart, double *tstop, double *dtout, int *algor, double *h0, double *tol, double jcen[], double *rcen, double *rmax, double en[], double am[], double *cefac, int *ndump, int *nfun, int *nbod, int *nbig, double m[], double x[][3], double v[][3], double s[][3], double rho[], double rceh[], int stat[], char id[][8], double ngf[][4], int opt[], int *opflag, int *ngflag, char outfile[][80], char dumpfile[][80], char mem[][80], int lmem[], void (*mdt_hy_) FUNCSIG_MAL_HVAR);


//Function for output
extern void mio_dump_(double *time, double *tstart, double *tstop, double *dtout, int *algor, double *h0, double *tol, double jcen[], double *rcen, double *rmax, double en[], double am[], double *cefac, int *ndump, int *nfun, int *nbod, int *nbig, double m[], double x[][3], double v[][3], double s[][3], double rho[], double rceh[], int stat[], char id[][8], double ngf[][4], double epoch[], int opt[], int *opflag, char dumpfile[][80], char mem[][80], int lmem[]);


//Other functions
extern void mxx_en_(double jcen[], int* nbod, int* nbig, double m[], double xh[][3], double vh[][3], double s[][3], double* e, double* l2);

extern void mxx_sync_(double* time, double* tstart, double* h0, double* tol, double jcen[], int* nbod, int* nbig, double m[], double x[][3], double v[][3], double s[][3], double rho[], double rceh[], int stat[], char id[][8], double epoch[], double ngf[][4], int opt[], int *ngflag);


int main()
{
	int j,algor,nbod,nbig, stat[NMAX], lmem[NMESS];
	int opflag,ngflag,ndump,nfun;
	double m[NMAX],xh[NMAX][3],vh[NMAX][3],s[NMAX][3],rho[NMAX];
	double rceh[NMAX],epoch[NMAX],ngf[NMAX][4],rmax,rcen,jcen[3];
	double cefac,time,tstart,tstop,dtout,h0,tol,en[3],am[3];
	char id[NMAX][8];
	char outfile[3][80], dumpfile[4][80], mem[NMESS][80];

	int opt[8] = {0,1,1,2,0,1,0,0};


	// Get initial conditions and integration parameters
	mio_in_(&time, &tstart, &tstop, &dtout, &algor, &h0, &tol, &rmax, &rcen, jcen, en, am, &cefac, &ndump, &nfun, &nbod, &nbig, m, xh, vh, s, rho, rceh, stat, id, epoch, ngf, opt, &opflag, &ngflag, outfile, dumpfile, lmem, mem);


	// If this is a new integration, integrate all the objects to a common epoch.
	if(opflag==-2)
	{
		//Ignoring some output
		mxx_sync_(&time, &tstart, &h0, &tol, jcen, &nbod, &nbig, m, xh, vh, s, rho, rceh, stat, id, epoch, ngf, opt, &ngflag);
		opflag = -1;
	}


	//Main Integration
	if(algor == 1) mal_hcon_(&time, &tstart, &tstop, &dtout, &algor, &h0, &tol, jcen, &rcen, &rmax, en, am, &cefac, &ndump, &nfun, &nbod, &nbig, m, xh, vh, s, rho, rceh, stat, id, ngf, opt, &opflag, &ngflag, outfile, dumpfile, mem, lmem, &mdt_mvs_, &mco_h2mvs_, &mco_mvs2h_);

	if(algor == 9) mal_hcon_(&time, &tstart, &tstop, &dtout, &algor, &h0, &tol, jcen, &rcen, &rmax, en, am, &cefac, &ndump, &nfun, &nbod, &nbig, m, xh, vh, s, rho, rceh, stat, id, ngf, opt, &opflag, &ngflag, outfile, dumpfile, mem, lmem, &mdt_mvs_, &mco_iden_, &mco_iden_);

	if(algor == 2) mal_hvar_(&time, &tstart, &tstop, &dtout, &algor, &h0, &tol, jcen, &rcen, &rmax, en, am, &cefac, &ndump, &nfun, &nbod, &nbig, m, xh, vh, s, rho, rceh, stat, id, ngf, opt, &opflag, &ngflag, outfile, dumpfile, mem, lmem, &mdt_bs1_);

	if(algor == 3) mal_hvar_(&time, &tstart, &tstop, &dtout, &algor, &h0, &tol, jcen, &rcen, &rmax, en, am, &cefac, &ndump, &nfun, &nbod, &nbig, m, xh, vh, s, rho, rceh, stat, id, ngf, opt, &opflag, &ngflag, outfile, dumpfile, mem, lmem, &mdt_bs2_);

	if(algor == 4) mal_hvar_(&time, &tstart, &tstop, &dtout, &algor, &h0, &tol, jcen, &rcen, &rmax, en, am, &cefac, &ndump, &nfun, &nbod, &nbig, m, xh, vh, s, rho, rceh, stat, id, ngf, opt, &opflag, &ngflag, outfile, dumpfile, mem, lmem, &mdt_ra15_);

	if(algor == 10) mal_hcon_(&time, &tstart, &tstop, &dtout, &algor, &h0, &tol, jcen, &rcen, &rmax, en, am, &cefac, &ndump, &nfun, &nbod, &nbig, m, xh, vh, s, rho, rceh, stat, id, ngf, opt, &opflag, &ngflag, outfile, dumpfile, mem, lmem, &mdt_hy_, &mco_h2dh_, &mco_dh2h_);

	if(algor == 11) mal_hcon_(&time, &tstart, &tstop, &dtout, &algor, &h0, &tol, jcen, &rcen, &rmax, en, am, &cefac, &ndump, &nfun, &nbod, &nbig, m, xh, vh, s, rho, rceh, stat, id, ngf, opt, &opflag, &ngflag, outfile, dumpfile, mem, lmem, &mdt_cb_, &mco_h2cb_, &mco_cb2h_);

	if(algor == 12) mal_hcon_(&time, &tstart, &tstop, &dtout, &algor, &h0, &tol, jcen, &rcen, &rmax, en, am, &cefac, &ndump, &nfun, &nbod, &nbig, m, xh, vh, s, rho, rceh, stat, id, ngf, opt, &opflag, &ngflag, outfile, dumpfile, mem, lmem, &mdt_wb_, &mco_h2wb_, &mco_wb2h_);


	// Do a final data dump
	for(j=2; j<=nbod; j++)
		epoch[j] = time;

	mio_dump_(&time, &tstart, &tstop, &dtout, &algor, &h0, &tol, jcen, &rcen, &rmax, en, am, &cefac, &ndump, &nfun, &nbod, &nbig, m, xh, vh, s, rho, rceh, stat, id, ngf, epoch, opt, &opflag, dumpfile, mem, lmem);


	// Calculate and record the overall change in energy and ang. momentum

	//Ignoring some output
	// ...

	mxx_en_(jcen, &nbod, &nbig, m, xh, vh, s, &(en[2]), &(am[2]));

	//Ignoring some output
	// ...

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
