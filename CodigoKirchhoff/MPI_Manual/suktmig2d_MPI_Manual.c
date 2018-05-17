#include "su.h"
#include "segy.h"
#include "header.h"
#include "avx_mathfun.h"
#include <signal.h>
#include <mpi.h>
#include <time.h>
#include <sys/time.h>
#include <immintrin.h>

char *sdoc[] = {
" 									",
" SUKTMIG2D - prestack time migration of a common-offset section with	",
"		the double-square root (DSR) operator			",
" 									",
"  suktmig2d < infile vfile= [parameters]  > outfile			",
"  									",
" 									",
" Required Parameters:							",
" vfile=	rms velocity file (units/s) v(t,x) as a function of time",
" dx=		distance (units) between consecutive traces		",
" nvelcdp       number of consecutive cdps in velocity file 	        ",
" firstcdp      first cdp number in velocity file			",
" lastcdp  	last cdp number in velocity file			",
" offmin=	minimum offset						",
" offmax=	maximum offset						",
"									",
" Optional parameters:							",
" intoff= 	interval between offsets				",
" cdp_trace_last                last cdp in data            		",
" angmax=40	maximum aperture angle for migration (degrees)		",
" hoffset=.5*tr.offset		half offset (m)				",
" nfc=16	number of Fourier-coefficients to approximate low-pass	",
"		filters. The larger nfc the narrower the filter		",
" fwidth=5 	high-end frequency increment for the low-pass filters	",
" 		in Hz. The lower this number the more the number of	",
"		lowpass filters to be calculated for each input trace.	",
"									",
" Notes:								",
" Data must be preprocessed with sufrac to correct for the wave-shaping	",
" factor using phasefac=.25 for 2D migration.				",
"									",
" Input traces must be sorted into offset and cdp number. The velocity	",
" file consists of rms velocities for all CMPs as a function of vertical",
" time and horizontal position v(t,z) in C-style binary floating point	",
" numbers. It's easiest to supply v(t,z) that has the same dimensions as",
" the input data to be migrated.					",
"									",
" The units may be feet or meters, as long as these are consistent for	",
" Antialias filter is performed using (Gray,1992, Geoph. Prosp), using	",
" nc low- pass filtered copies of the data. The cutoff frequencies are	",
" calculated  as fractions of the Nyquist frequency.			",
"									",
" The maximum allowed angle is 80 degrees(a 10 degree taper is applied	",
" to the end of the aperture)						",
NULL};


/*
 * Credits:
 * Serial Version: CWP, Baoniu Han, bhan@dix.mines.edu, April 19th, 1998
 * MPI Version: CPGG/UFBA, Reynam Pestana, reynam@cpgg.ufba.br
 *
 * Trace header fields accessed: ns, dt, delrt, d2
 * Trace header fields modified: ns, dt, delrt
 */

/**************** end self doc *******************************************/

/* Prototypes for functions used internally */
#define LOOKFAC 2       /* Look ahead factor for npfaro   */
#define PFA_MAX 720720  /* Largest allowed nfft           */

// #define TAG_REQUEST  		10
// #define TAG_RESPONSE 		20
#define TAG_WORK_REQUEST	10
// #define TAG_WORK_REPLY		20
#define TAG_SEND_WORK		30
#define TAG_SEND_WORK_INTEGERS	31
#define TAG_SEND_WORK_FLOATS	32
#define TAG_SEND_WORK_DATA	33
#define TAG_REPLY		20
#define TAG_REPLY_INTEGERS	21
#define TAG_REPLY_DATA		22
#define TAG_FINISHED 		99
#define TAG_PARAMETERS		80

#define MAX_ELEMENTS		1024

MPI_Status status;
MPI_Request request;

FILE *in_file;
FILE *out_file;

/* Prototype of functions used internally */
void lpfilt(int nfc, int nfft, float dt, float fhi, float *filter);
void lpfilt_vec(int nfc, int nfft, float dt, float fhi, float *filter);
void migkt2d(float **data, int ntr, int nt, float dx, float dt, float tmax,
             int nfft, float fnyq, float h, float *fc, int nf, int nc, int nfc,
             int cdp_trace_first, float angmax, float **vel, float **mig, int dxcdp, 
             int firstcdp);
void gravar_dados(float **mig, int ntr, int nt, int oldoffset, char *arq);
struct timeval inicia_tempo();
double para_tempo(struct timeval tv1);
double tempo_offset;
void gravar_tempo_offset(int oldoffset, double tempo, char *arq);

segy intrace; 	/* input traces */
segy outtrace;	/* migrated output traces */

void gravar_velocidade(float **vel, int nvelcdp, int nt, int oldoffset, char *arq)
{
	FILE *arquivo;
	warn("Offset=%d, nvelcdp=%d, nt=%d", oldoffset, nvelcdp, nt);
	arquivo = fopen(arq,"a");
	int ix, it;
	for (ix=0; ix<nvelcdp; ++ix) 
	{
		fprintf(arquivo,"Inicio ix=%d - Offset %d\n",ix,oldoffset);
		for (it=0; it<nt; ++it) {
			fprintf(arquivo,"%f ",vel[ix][it]);
		}
		fprintf(arquivo,"\nFim ix=%d - Offset %d\n",ix,oldoffset);
		
	}
	fclose(arquivo);	
}

void get_parameters(segy *trace, int *nt, float *dt, float *tmax,  
	int *nfft, int *nf, float *fnyq, int *nc, float ***vel, float **fc, int fwidth, 
	int nvelcdp, int verbose, char *vfile)
{
	int i;

	FILE *vfp=NULL;

	*nt=intrace.ns;
	*dt=(float)intrace.dt/1000000;
	*tmax=(*nt-1)*(*dt);
	// *gottrace=1;

	/* Set up FFT parameters */
	*nfft = npfaro(*nt, LOOKFAC*(*nt));

	if(*nfft>= SU_NFLTS || *nfft >= PFA_MAX)
	  err("Padded nt=%d -- too big",*nfft);

	*nf = (*nfft)/2 + 1;

	/* Determine number of filters for antialiasing */
	*fnyq= 1.0/(2*(*dt));
	*nc=ceil(*fnyq/fwidth);
	if (verbose)
		warn(" The number of filters for antialiasing is nc= %d",*nc);

	*vel =   alloc2float(*nt,nvelcdp);
	*fc  =	alloc1float(*nc+1);

	for (i=1; i<*nc+1; ++i) {
		(*fc)[i]= *fnyq*i/(*nc);
	}

	vfp=efopen(vfile,"r");
	efread((*vel)[0],FSIZE,(*nt)*nvelcdp,vfp);
	efclose(vfp);

	warn("Vel file readed");
}

int read_offset(FILE *hfp, int *hfp_index, FILE *list_hfp[], int *ntr, 
	int offmin, int offmax, int *noffset, int *lastoffset,
	int intoff, int *cdp_trace_first, float *h, int *deltacdp, int nt, 
	int verbose, int nvelcdp, float ***data, int *oldoffset)
{
	int gottrace=1;
	int i,it,ix;	/* counters */
	// int lastoffset = 0;
//	int oldoffset = 0, 
int curoffset = 0;
	int olddeltacdp = 0;
	int oldcdp = 0;	/* temporary storage			*/
	int cdp_trace_last = 0;	/* last cdp in data		*/
	int ncdp = 0;	/* number of cdps in the velocity file	*/

	cwp_Bool check_cdp=cwp_false;	/* check cdp in velocity file	*/

	FILE *tracefp = NULL;	/* temp file to hold traces*/

	// While has more traces to migration from input file
	// for(;gottrace==1;){
		tracefp = etmpfile();
		hfp = etmpfile();
		*ntr = 0;

		curoffset=intrace.offset;
		do {
			*oldoffset = curoffset;
			curoffset = intrace.offset;
			if(*oldoffset != curoffset){
				if(curoffset>offmax){
					gottrace=2;
					break;
				}
				if((*noffset != 0) && (*ntr > 0))
					break;
			}
			if (curoffset >= offmin){
				if(*ntr == 0){
					if(*noffset > 0){
						if(curoffset != (*lastoffset + intoff)){
							if(!fgettr(in_file, &intrace)){
								gottrace=2;
								break;
							}
							continue;
						}
						else{
							*lastoffset = intrace.offset;
							(*noffset)++;
							*cdp_trace_first = intrace.cdp;
							*h=.5*intrace.offset;
						}

					}
					else{
						*lastoffset = intrace.offset;
						(*noffset)++;
						*cdp_trace_first = intrace.cdp;
						*h=.5*intrace.offset;
					}
				}
				++(*ntr);
				/* get new deltacdp value */
				*deltacdp=intrace.cdp-oldcdp;
				/* read headers and data */
				efwrite(&intrace,HDRBYTES, 1, hfp);
				efwrite(intrace.data, FSIZE, nt, tracefp);

				/* error trappings. */
				/* ...did cdp value interval change? */
				if ((*ntr>3) && (olddeltacdp!=*deltacdp)) {
					if (verbose) {
						warn("cdp interval changed in data");
						warn("ntr=%d olddeltacdp=%d deltacdp=%d",*ntr,olddeltacdp,*deltacdp);
					 	check_cdp=cwp_true;
					}
				}

				/* save cdp and deltacdp values */
				oldcdp=intrace.cdp;
				olddeltacdp=*deltacdp;

				cdp_trace_last = oldcdp;

			}
			if(!fgettr(in_file, &intrace))
				gottrace=2;
		}
		while (gottrace!=2);
		//warn("Saiu do gottrace!=2");
		if(*noffset==1){
			if (intoff==0)
				intoff = curoffset - *oldoffset;
			else
				if (((curoffset - *oldoffset)>intoff)||((intoff%(curoffset-*oldoffset))!=0)){
					gottrace = 0;
					warn("Invalid interval");
					return gottrace;
					// break;
				}
			}

		if(*ntr==0){
			gottrace = 0;
			// warn ("Finishing Process");
			return gottrace;
			// break;
		}

		ncdp=cdp_trace_last-*cdp_trace_first+1;

		/* error trappings */
		if ( (*cdp_trace_first==cdp_trace_last)
			|| (*deltacdp==0)
			|| (check_cdp==cwp_true) ) warn("Check cdp values in data!");

		erewind(tracefp);
		erewind(hfp);

		if ( (ncdp > nvelcdp) )
			warn("Check ncdp values in data great than ncdpvel ");

		//warn("0:alocando %d * %d", nt, ntr);
		*data = 	alloc2float(nt,*ntr);
		
		for (ix=0; ix<*ntr; ++ix)
			efread((*data)[ix],FSIZE,nt,tracefp);

		efclose(tracefp);
		//mig = alloc2float(nt,ntr);

		warn("Offset %d",*oldoffset);

		/* Sending work */
		// hfp_index = work_sent % MAX_ELEMENTS;
		list_hfp[*hfp_index] = hfp;

		// MPI_Recv(&msg, 0, MPI_INT, MPI_ANY_SOURCE, 
		// 	TAG_WORK_REQUEST, MPI_COMM_WORLD, &status);
		// origin = status.MPI_SOURCE;

		// warn("0:Work request from [%d]", origin);

		//warn("[ntr %d, nt %d]", ntr, nt);
		// integers_data[0] = ntr;				// 0. int ntr=0;
		// integers_data[1] = nt;				// 1. int nt;
		// integers_data[2] = nfft;			// 2. int nfft;
		// integers_data[3] = nf;				// 3. int nf;
		// integers_data[4] = nc;				// 4. int nc;
		// integers_data[5] = nfc;				// 5. int nfc;
		// integers_data[6] = cdp_trace_first;		// 6. int cdp_trace_first=0;
		// integers_data[7] = deltacdp;			// 7. int deltacdp=0
		// integers_data[8] = firstcdp;			// 8. int firstcdp=0;
		// integers_data[9] = hfp_index;
		// warn("0:Enviando dados inteiros");
		// MPI_Send(&integers_data, 9, MPI_INT, 
			// origin, TAG_REQUEST, MPI_COMM_WORLD);

		// floats_data[0] = dx;				//  9. float dx=0.0;
		// floats_data[1] = dt;				// 10. float dt;
		// floats_data[2] = tmax;				// 11. float tmax;
		// floats_data[3] = fnyq;				// 12. float fnyq;
		// floats_data[4] = h;				// 13. float h;
		// floats_data[5] = angmax;			// 14. float angmax;
		// warn("0:Enviando dados de pontos flutuantes");
		// MPI_Send(&floats_data, 6, MPI_FLOAT, 
			// origin, TAG_REQUEST, MPI_COMM_WORLD);
		
								// 15. float *fc=NULL;
		// warn("0:Enviando dados de um vetor");
		// MPI_Send(fc, nc + 1, MPI_FLOAT, 
			// origin, TAG_REQUEST, MPI_COMM_WORLD);
		
								// 17. float **data=NULL;
		// warn("0:Enviando dados de uma matriz com os dados");
		// MPI_Send(data[0], ntr * nt, MPI_FLOAT, 
			// origin, TAG_REQUEST, MPI_COMM_WORLD);
		
								// 18. float **vel=NULL;
		// warn("0:Enviando dados de uma matriz de velocidades");
		// MPI_Send(vel[0], nvelcdp * nt, MPI_FLOAT, 
			// origin, TAG_REQUEST, MPI_COMM_WORLD);
		//MPI_Buffer_detach(i_old_buffer, &i_old_size);
								// 19. float **mig=NULL;
		/*MPI_Recv(mig[0], ntr * nt, MPI_FLOAT, (work_sent % slaves) + 1, 
			MPI_ANY_TAG, MPI_COMM_WORLD, &status);*/
		
		// warn("[work_sent %d]", work_sent);
		// work_sent++;
		// if (gottrace != 2)
		// {
			// do
			// {
				// if (flag)
				// {
					// warn("0:Preparing to receive response");
					//warn("IRecv [work_received: %d] ... ", work_received);
					/*
					mig = alloc2float(nt, ntr);
					MPI_Irecv(mig[0], ntr * nt, MPI_FLOAT, (work_received % slaves) + 1, 
						TAG_RESPONSE, MPI_COMM_WORLD, &request);
					*/
					// MPI_Irecv(&integers_data_return, 3, MPI_INT, MPI_ANY_SOURCE, 
						// TAG_PREPARE_RESPONSE, MPI_COMM_WORLD, &request);
					
					// flag = 0;
				// }
				// warn("0:Before if there is workers available");
				// if (work_sent - work_received >= slaves) MPI_Wait(&request, &status);
				// warn("0:After if there is workers available");
				// MPI_Test(&request, &flag, &status);
				// if (flag)
				// {
					// warn("0:Getting a response");
					// origin 			= status.MPI_SOURCE;
					// tag 			= status.MPI_TAG;
					// warn("0:origin=%d tag=%d", origin, tag);
					// ntr			= integers_data[0];
					// nt			= integers_data[1];
					// hfp_index		= integers_data[2];
					
					// mig = alloc2float(nt, ntr);
					// MPI_Recv(mig[0], ntr * nt, MPI_FLOAT, origin, 
						// TAG_RESPONSE, MPI_COMM_WORLD, &status);

					//warn("Processing [work_received: %d] ... ", work_received);
					//FILE *inner_hfp = list_hfp[work_received % MAX_ELEMENTS];
					// FILE *inner_hfp = list_hfp[hfp_index];
					// for (ix = 0; ix < ntr; ++ix)
					// {
						// efread(&outtrace, HDRBYTES, 1, inner_hfp);
						// for (it = 0; it < nt; ++it) 
						// {
							// outtrace.data[it] = mig[ix][it];
						// }
						// fputtr(out_file, &outtrace);
					// }
					// efclose(inner_hfp);
					//gravar_dados(mig,ntr,nt,offmin+(20*work_received),"dados_migrados_mpi.txt");
					// free2float(mig);
					// work_received++;
					// flag = -1;
				// }
			// }
			// while (flag);
		// }
	
		/*}*/
		
		/*migkt2d(data, ntr, nt, dx, dt, tmax, nfft, fnyq, h, fc, nf, nc, nfc,
      		        cdp_trace_first, angmax, vel, mig, deltacdp, firstcdp);*/

		/*
		for (ix=0; ix<ntr; ++ix) {
			efread(&outtrace, HDRBYTES, 1, hfp);
			for (it=0; it<nt; ++it) {
				outtrace.data[it] = mig[ix][it];
			}
			puttr(&outtrace);
		}
		efclose(hfp);
		//free2float(mig);
		*/
		// free2float(data);

		// if (gottrace==2){
		// 	warn("Finishing Process");
		// 	gottrace=0;
		// 	break;
		// }
	// } // for(;*gottrace==1;){
	return gottrace;
}

int main(int argc,char **argv){
	int size, rank, slaves, work_sent = 0, work_received = 0;
	//int rows, control;

	// int gottrace=1;
	int i,it,ix;	/* counters */
	int ntr=0,nt;	/* x,t */

	int verbose=0;	/* is verbose?				*/
	int nc;		/* number of low-pass filtered versions	*/
			/*  of the data for antialiasing	*/
	int nfft,nf;	/* number of frequencies		*/
	int nfc;	/* number of Fourier coefficients for low-pass filter */
	int fwidth;	/* high-end frequency increment for the low-pass */
			/* filters 				*/
	int firstcdp=0;	/* first cdp in velocity file		*/
	int lastcdp=0;	/* last cdp in velocity file		*/

	int cdp_trace_first=0;	/* first cdp in data		*/
	// int cdp_trace_last=0;	/* last cdp in data		*/
	int nvelcdp;     /* number of cdps in the velocity file */

	int offmin,offmax;
	int intoff;
	int lastoffset=0;

	int oldoffset=0;
//, curoffset=0;
	int noffset=0;
	// int oldcdp=0;	/* temporary storage			*/
	// int olddeltacdp=0;
	int deltacdp=0;
	// int ncdp=0;	/* number of cdps in the velocity file	*/

	float dx=0.0;	/* cdp sample interval */
	float hoffset=0.0;  /* half receiver-source */
			/* no aliasing of the operator */
	float dt;	/* t sample interval */
	float h;	/* offset */

	float angmax;   /* maximum aperture angle */

	float tmax;	/* maximum time */

	float fnyq;	/* Nyquist frequency */

	float *fc=NULL;		/* cut-frequencies for low-pass filters */

	float **vel=NULL;	/* array of velocity values from vfile */
	float **data=NULL;	/* input data array*/
	float **mig=NULL;	/* output migrated data array */

	/* file names */
	char *vfile="";		/* name of velocity file */
	char *input_file="";
	char *output_file="";
	// FILE *vfp=NULL;
	//FILE *tracefp=NULL;	/* temp file to hold traces*/
	FILE *hfp=NULL;		/* temp file to hold trace headers */

	cwp_Bool check_cdp=cwp_false;	/* check cdp in velocity file	*/

	FILE *list_hfp[MAX_ELEMENTS];
	int flag = -1;
	int integers_data[5];
	int integers_data_return[2];
	float floats_data[1];
	int msg = 0;
	int origin = 0 ;
	int tag = 0;
	int hfp_index = 0;
	int worker_finished = 0;

    	/*Inicialização do MPI*/
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	slaves = size-1;

	initargs(argc, argv);
	//requestdoc(1);
	requestdoc(0);

	MUSTGETPARFLOAT("dx",&dx);
	MUSTGETPARINT("nvelcdp",&nvelcdp);
	MUSTGETPARINT("firstcdp",&firstcdp);
	MUSTGETPARINT("lastcdp",&lastcdp);
	MUSTGETPARINT("offmin",&offmin);
	MUSTGETPARINT("offmax",&offmax);
	MUSTGETPARSTRING("vfile",&vfile);
	MUSTGETPARSTRING("input_file",&input_file);
	MUSTGETPARSTRING("output_file",&output_file);
	
	// warn("%d:nvelcdp %d", rank, nvelcdp);

	if (!getparfloat("angmax",&angmax)) angmax=40;
	if (!getparint("nfc",&nfc)) nfc=16;
	if (!getparint("fwidth",&fwidth)) fwidth=5;
	if (!getparint("verbose",&verbose)) verbose=0;
	if (!getparint("intoff",&intoff)) intoff=0;
	//warn("Antes do if rank==0");

	/* Processo Mestre */
	if (rank == 0) 
	{
		int has_more_offset;

		in_file = fopen(input_file, "r");
		// in_file = fopen("/home/rprado/dev/kirchhoff/MIG-TEMPO-OK/SYNT-DATA/synt-off.su" , "r");
		// in_file = fopen("/home/rodrigo/ic/dev/MIG-TEMPO-OK/REAL-DATA/jeqoffset_new.su" , "r");
		
		out_file = fopen(output_file, "w");
		// out_file = fopen("saida-synt-mpi.su" , "w");
		//FILE *arquivo = fopen("dados_migrados_mpi.txt", "w");
		//fclose(arquivo);
		
		if (!fgettr(in_file, &intrace)) err("can't get first trace");

		/* Sending header of first trace to workers */
		for (i = 0; i < slaves; i++) {
			MPI_Send(&intrace, HDRBYTES, MPI_BYTE, i+1, TAG_PARAMETERS, MPI_COMM_WORLD);
		}

		get_parameters(&intrace, &nt, &dt, &tmax, &nfft, &nf, &fnyq, &nc, 
			&vel, &fc, fwidth, nvelcdp, verbose, vfile);

		/*warn("0:nt=%d dt= %f tmax=%f nfft=%d nf=%d fnyq=%f nc=%d", 
			nt, dt, tmax, nfft, nf, fnyq, nc);*/

		// gravar_velocidade(vel,nvelcdp,nt,oldoffset,"vel_file_gravado.txt");

		has_more_offset = read_offset(hfp, &hfp_index, list_hfp, &ntr, 
			offmin, offmax, &noffset, &lastoffset, intoff,
			&cdp_trace_first, &h, &deltacdp, nt, verbose, nvelcdp, &data, &oldoffset);
	
		/*warn("0:hfp_index=%d ntr= %d noffset=%d lastoffset=%d cdp_trace_first=%d h=%f deltacdp=%d", 
			hfp_index, ntr, noffset, lastoffset, cdp_trace_first, h, deltacdp);*/

		/* Handle messages */

		/*while (has_more_offset != 0 && (work_sent > work_received 
			|| worker_finished < slaves))*/
		while (has_more_offset != 0 || worker_finished < slaves)
		{
			MPI_Recv(&msg, 0, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, 
				MPI_COMM_WORLD, &status);

			origin = status.MPI_SOURCE;
			tag = status.MPI_TAG;

			warn("0:Processing message from [%d] tag [%d]", origin, tag);

			if (tag == TAG_WORK_REQUEST)
			{
				if (has_more_offset == 0)
				{
					warn("0:Sending a message to finish to [%d]", origin, tag);
					MPI_Send(&msg, 0, MPI_INT, origin, 
						TAG_FINISHED, MPI_COMM_WORLD);
					++worker_finished;
				}
				else // gottrace == 1 and gottrace == 2 
				{
					MPI_Send(&msg, 0, MPI_INT, origin, 
						TAG_SEND_WORK, MPI_COMM_WORLD);

					// warn("0:Enviando dados inteiros");
					integers_data[0] = ntr;				// 0. int ntr=0;
					integers_data[1] = cdp_trace_first;		// 6. int cdp_trace_first=0;
					integers_data[2] = deltacdp;			// 7. int deltacdp=0
					integers_data[3] = hfp_index;
integers_data[4] = oldoffset;
					MPI_Send(&integers_data, 5, MPI_INT, 
						origin, TAG_SEND_WORK_INTEGERS, 
						MPI_COMM_WORLD);
					
					// integers_data[1] = nt;				// 1. int nt;
					// integers_data[2] = nFimft;			// 2. int nfft;
					// integers_data[3] = nf;				// 3. int nf;
					// integers_data[4] = nc;				// 4. int nc;
					// integers_data[5] = nfc;				// 5. int nfc;
					// integers_data[8] = firstcdp;			// 8. int firstcdp=0;
					
					// warn("0:Enviando dados de pontos flutuantes");
					floats_data[0] = h;				// 13. float h;
					MPI_Send(&floats_data, 1, MPI_FLOAT, 
						origin, TAG_SEND_WORK_FLOATS, 
						MPI_COMM_WORLD);
					
					// floats_data[0] = dx;				//  9. float dx=0.0;
					// floats_data[1] = dt;				// 10. float dt;
					// floats_data[2] = tmax;				// 11. float tmax;
					// floats_data[3] = fnyq;				// 12. float fnyq;
					// floats_data[5] = angmax;			// 14. float angmax;
					
											// 17. float **data=NULL;
					// warn("0:Enviando dados de uma matriz com os dados");
					MPI_Send(data[0], ntr * nt, MPI_FLOAT, 
						origin, TAG_SEND_WORK_DATA, 
						MPI_COMM_WORLD);
					

					++work_sent;

					if (has_more_offset == 2) has_more_offset = 0;
					else
					{
						/* Read a new offset to send */

						hfp_index++;

						has_more_offset = read_offset(hfp, &hfp_index, list_hfp, &ntr, 
							offmin, offmax, &noffset, &lastoffset, intoff,
							&cdp_trace_first, &h, &deltacdp, nt, verbose, nvelcdp, &data, &oldoffset);
						/*warn("0:has_more_offset=%d", has_more_offset);
						warn("0:hfp_index=%d ntr= %d oldoffset=%d curoffset=%d noffset=%d  lastoffset=%d cdp_trace_first=%d h=%f deltacdp=%d", 
							hfp_index, ntr, noffset, lastoffset, cdp_trace_first, h, deltacdp);*/

					}
				}
			}
			else if (tag == TAG_REPLY)
			{
				int inner_ntr;
				int inner_hfp_index;
				/* Receive migrated data */

				warn("0:Receiving integers");

				MPI_Recv(&integers_data_return, 2, MPI_INT, origin, 
					TAG_REPLY_INTEGERS, MPI_COMM_WORLD, &status);
				inner_ntr		= integers_data_return[0];
				// nt			= integers_data[1];
				inner_hfp_index		= integers_data_return[1];
				warn("0:inner_hfp_index=%d", inner_hfp_index);
				warn("0:Receiving data");

				mig = alloc2float(nt, inner_ntr);
				MPI_Recv(mig[0], inner_ntr * nt, MPI_FLOAT, origin, 
					TAG_REPLY_DATA, MPI_COMM_WORLD, &status);
				
				//FILE *inner_hfp = list_hfp[work_received % MAX_ELEMENTS];
				FILE *inner_hfp = list_hfp[inner_hfp_index];
				for (ix = 0; ix < inner_ntr; ++ix) 
				{
					efread(&outtrace, HDRBYTES, 1, inner_hfp);
					for (it = 0; it < nt; ++it) 
					{
						outtrace.data[it] = mig[ix][it];
					}
					fputtr(out_file, &outtrace);
				}
				efclose(inner_hfp);
				// gravar_dados(mig,ntr,nt,offmin+(20*work_received),"dados_migrados_mpi.txt");
				free2float(mig);
				++work_received;
			}

		} // while (has_more_offset != 0 && (work_sent > work_received 

		// if (!flag)
		// {
		// 	//warn("Waiting [work_received: %d] ... ", work_received);
		// 	MPI_Wait(&request, &status);

		// 	origin 			= status.MPI_SOURCE;

		// 	ntr			= integers_data[0];
		// 	nt			= integers_data[1];
		// 	hfp_index		= integers_data[2];
			
		// 	mig = alloc2float(nt, ntr);
		// 	MPI_Recv(mig[0], ntr * nt, MPI_FLOAT, origin, 
		// 		TAG_RESPONSE, MPI_COMM_WORLD, &status);


		// 	// war\n("Processando trabalho pendente [%d] ... ", work_received);

		// 	// mig = alloc2float(nt,ntr);
		// 	/*
		// 	MPI_Recv(mig[0], ntr * nt, MPI_FLOAT, (work_received % slaves) + 1, 
		// 		TAG_RESPONSE, MPI_COMM_WORLD, &status);
		// 	*/
		// 	//FILE *inner_hfp = list_hfp[work_received % MAX_ELEMENTS];
		// 	FILE *inner_hfp = list_hfp[hfp_index];
		// 	for (ix = 0; ix < ntr; ++ix) 
		// 	{
		// 		efread(&outtrace, HDRBYTES, 1, inner_hfp);
		// 		for (it = 0; it < nt; ++it) 
		// 		{
		// 			outtrace.data[it] = mig[ix][it];
		// 		}
		// 		fputtr(out_file, &outtrace);
		// 	}
		// 	efclose(inner_hfp);
		// 	// gravar_dados(mig,ntr,nt,offmin+(20*work_received),"dados_migrados_mpi.txt");
  //   			free2float(mig);
  //   			work_received++;
  //   			//warn("[%d] trabalhos pendente processado.", work_received);
		// }

		// for (; work_received < work_sent; work_received++)
		// {
		// 	//warn("Processando trabalho finais [%d] ... ", work_received);
		// 	//warn("Recv [work_received: %d] ... ", work_received);

		// 	MPI_Recv(&integers_data_return, ntr * nt, MPI_FLOAT, origin, 
		// 		TAG_PREPARE_RESPONSE, MPI_COMM_WORLD, &status);
			
		// 	origin 			= status.MPI_SOURCE;

		// 	ntr			= integers_data[0];
		// 	nt			= integers_data[1];
		// 	hfp_index		= integers_data[2];
			
		// 	mig = alloc2float(nt, ntr);
		// 	MPI_Recv(mig[0], ntr * nt, MPI_FLOAT, origin, 
		// 		TAG_RESPONSE, MPI_COMM_WORLD, &status);
		// 	/*
		// 	mig = alloc2float(nt, ntr);
		// 	MPI_Recv(mig[0], ntr * nt, MPI_FLOAT, (work_received % slaves) + 1, 
		// 		TAG_RESPONSE, MPI_COMM_WORLD, &status);
		// 	*/
		// 	//FILE *inner_hfp = list_hfp[work_received % MAX_ELEMENTS];
		// 	FILE *inner_hfp = list_hfp[hfp_index];
		// 	for (ix = 0; ix < ntr; ++ix) 
		// 	{
		// 		efread(&outtrace, HDRBYTES, 1, inner_hfp);
		// 		for (it = 0; it < nt; ++it) 
		// 		{
		// 			outtrace.data[it] = mig[ix][it];
		// 		}
		// 		fputtr(out_file, &outtrace);
		// 	}

		// 	efclose(inner_hfp);
		// 	// gravar_dados(mig,ntr,nt,offmin+(20*work_received),"dados_migrados_mpi.txt");
		// 	free2float(mig);
		// 	//work_received++;

		// 	//warn("[%d] trabalhos processados (final).", work_received+1);
		// }

		// //warn("Enviando mensagens de termino.");

		// /*for (i = 0; i < slaves; i++)
		// 	MPI_Send(&integers_data, 0, MPI_INT, i+1, TAG_FINISHED, MPI_COMM_WORLD);*/
		// for (i = 0; i < slaves; i++) {
		// 	MPI_Recv(&msg, 0, MPI_INT, i+1, TAG_WORK_REQUEST, MPI_COMM_WORLD, &status);
		// 	MPI_Send(&integers_data, 0, MPI_INT, i+1, TAG_FINISHED, MPI_COMM_WORLD);
		// }

		warn("0:Finishing Process");

		fclose(in_file);
		fclose(out_file);
	}
	else // Slaves
	{
		//warn("Rank %d", rank);
		// int cut_frequencies_data[9];
		// int input_data[9];
		// int input_velocities_data[9];
		// int output_data[9];
		
		// MPI_Request request;
		int tag = 0;
		//MPI_Status status;

		/* Receiving header from master */
		MPI_Recv(&intrace, HDRBYTES, MPI_BYTE, 0, MPI_ANY_TAG, 
			MPI_COMM_WORLD, &status);
    		
		get_parameters(&intrace, &nt, &dt, &tmax, &nfft, &nf, &fnyq, &nc, 
			&vel, &fc, fwidth, nvelcdp, verbose, vfile);
		
		/*warn("%d:nt=%d dt= %f tmax=%f nfft=%d nf=%d fnyq=%f nc=%d", 
			rank, nt, dt, tmax, nfft, nf, fnyq, nc);*/

		/* Requesting work */
		warn("%d:Resquesting work", rank);
		MPI_Send(&msg, 0, MPI_INT, 0, TAG_WORK_REQUEST, MPI_COMM_WORLD);
		
		/* Getting response */
		MPI_Recv(&msg, 0, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		tag = status.MPI_TAG;

		while (tag != TAG_FINISHED)
		{
			/* Receiving data */
			//warn("%d:Arguardando requisicao", rank);
			warn("%d:Receiving integers data", rank);
			MPI_Recv(&integers_data, 5, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			tag = status.MPI_TAG;
			
			ntr			= integers_data[0];
			cdp_trace_first		= integers_data[1];
			deltacdp		= integers_data[2];
			hfp_index		= integers_data[3];
oldoffset = integers_data[4];
				// nt			= integers_data[1];
				// nfft 			= integers_data[2];
				// nf			= integers_data[3];
				// nc			= integers_data[4];
				// nfc			= integers_data[5];
				// firstcdp		= integers_data[8];
			warn("%d:hfp_index=%d", rank, hfp_index);
			warn("%d:Receiving floats data", rank);
			MPI_Recv(&floats_data, 1, MPI_FLOAT, 0, MPI_ANY_TAG, 
				MPI_COMM_WORLD, &status);
			h			= floats_data[0];
				// dx			= floats_data[0];
				// dt			= floats_data[1];
				// tmax 			= floats_data[2];
				// fnyq 			= floats_data[3];
				// angmax 			= floats_data[5];

				//warn("%d:Recebendo dados um vetor", rank);
				// fc = alloc1float(nc + 1);
				// MPI_Recv(fc, nc + 1, MPI_FLOAT, 0, MPI_ANY_TAG, 
					// MPI_COMM_WORLD, &status);

			warn("%d:Receiving data", rank);
			data = alloc2float(nt, ntr);
			MPI_Recv(data[0], ntr * nt, MPI_FLOAT, 0, MPI_ANY_TAG, 
				MPI_COMM_WORLD, &status);
				
				//warn("%d:Recebendo dados de velocidades", rank);
				// vel = alloc2float(nt, nvelcdp);
				// MPI_Recv(vel[0], nvelcdp * nt, MPI_FLOAT, 0, MPI_ANY_TAG, 
					// MPI_COMM_WORLD, &status);

			/* Processing data */
			
//char nome_arq[255];
//tempo_offset = 0;
//sprintf(nome_arq, "mpi_offset_time_%d.txt", oldoffset);
			warn("%d:Processing work", rank);
			mig = alloc2float(nt, ntr);
			migkt2d(data, ntr, nt, dx, dt, tmax, nfft, fnyq, h, fc, nf, nc, nfc, 
				cdp_trace_first, angmax, vel, mig, deltacdp, firstcdp);
//gravar_tempo_offset(oldoffset,tempo_offset,nome_arq);
				
			/* Sending data */

			MPI_Send(&msg, 0, MPI_INT, 0, TAG_REPLY, MPI_COMM_WORLD);

			integers_data_return[0] = ntr;
			// integers_data_return[1] = nt;
			integers_data_return[1] = hfp_index;
			MPI_Send(&integers_data_return, 2, MPI_INT, 0, TAG_REPLY_INTEGERS, MPI_COMM_WORLD);

			//warn("%d:Enviando resposta", rank);
			MPI_Send(mig[0], ntr * nt, MPI_FLOAT, 0, TAG_REPLY_DATA, MPI_COMM_WORLD);

				/* For debug purposes */
				/*warn("%d:%d:%d:%d:%d:%d:%d:%d:%d:%d", rank, integers_data[0], 
					integers_data[1], integers_data[2], integers_data[3], 
					integers_data[4], integers_data[5], integers_data[6],
					integers_data[7], integers_data[8]);
				warn("%d:%f:%f:%f:%f:%f:%f", rank, floats_data[0], floats_data[1],
					floats_data[2], floats_data[3], floats_data[4], floats_data[5]);
				for (i = 0; i <= nc; i++) 
				{
					warn("%d:%f", rank, fc[i+1]);
				}
				
				for (i = 0; i < ntr; ++i)
					for (ix = 0; ix < nt; ++ix) 
						mig[i][ix] = 5.5; 
				
				for (i = 0; i < ntr; ++i)
					for (ix = 0; ix < nt; ++ix)
						warn("%d:%d:%d:%f:%f:%f", rank, i, ix, data[i][ix], 
							vel[i][ix], mig[i][ix]);
				*/

			// Freeing resources
			// free1float(fc);
			free2float(data);
			free2float(mig);
			// free2float(vel);
	        	

			// MPI_Recv(&control, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
			// MPI_Recv(&control, 2, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);

			//MPI_Wait(&request, &status);
			//MPI_Recv(&rows, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
			//MPI_Recv(&a, rows * N, MPI_LONG, 0, 1, MPI_COMM_WORLD, &status);
			//MPI_Irecv(&b, N * N, MPI_LONG, 0, 1, MPI_COMM_WORLD, &request);
			//MPI_Wait(&request, &status);

			/* Requesting work */
			warn("%d:Resquesting work", rank);
			MPI_Send(&msg, 0, MPI_INT, 0, TAG_WORK_REQUEST, MPI_COMM_WORLD);
			
			/* Getting response */
			MPI_Recv(&msg, 0, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			tag = status.MPI_TAG;
		
		} // while (tag != TAG_FINISHED)
	} // else // Slaves

	int retorno = CWP_Exit();
	MPI_Finalize(); //finalizar MPI
	return(retorno);
}

void  migkt2d ( float **data, int ntr, int nt, float dx, float dt, float tmax, 
               int nfft, float fnyq, float h, float *fc, int nf, int nc, int nfc, 
               int cdp_trace_first,float angmax, float **vel,float **mig, int dxcdp, int firstcdp)   

{
	int k,imp,iip,it,ifc;	/* counters */

	float p=0.0;	/* horizontal slowness of the migration operator */
	/*float pmin = 0.0;*/		/* maximum horizontal slowness for which there's */
			/* no aliasing of the operator */
	float x;	/* aperture distance */
	float xmax=0.0;	/* maximum aperture distance */

	float obliq;	/* obliquity factor */
	float geoms;	/* geometrical spreading factor */

	float mp,ip;	/* mid-point and image-point coordinates */
	float t;	/* time */
	float t0;	/* vertical traveltime */

	/*float ang;*/	/* aperture angle */
	float angtaper = 0.0;	/* aperture-angle taper */
	float v;		/* velocity */
  
	float *filter = NULL;	/* array of low-pass filter values */

	float **lowpass = NULL; /* low-pass filtered version of the trace */

	// register float *rtin=NULL,*rtout=NULL;/* real traces */
	// register complex *ct=NULL;   /* complex trace */
	float *rtin=NULL,*rtout=NULL;/* real traces */
	complex *ct=NULL;   /* complex trace */
	//float rtin[nfft] __attribute__((aligned(16)));;
	//float rtout[nfft] __attribute__((aligned(16)));;
	
	float datalo[8]; 		//__attribute__((aligned(16)));;
	float datahi[8]; 		//__attribute__((aligned(16)));;
	int itb, ite;
	float firstt, amplo, amphi;

	float *t_vec=NULL;
	float *tr_vec=NULL;
	float *ts_vec=NULL;
	float *geoms_vec=NULL;
	float *obliq_vec=NULL;
	float *pmin_vec=NULL;
	float p_vec[ntr+1];
	float fplo_vec[ntr+1];
	/*float fphi_vec[ntr+1];*/
	float ref_arg_vec[ntr+1];
	/*float wlo_vec[ntr+1];*/
	/*float whi_vec[ntr+1];*/
	float itb_vec[ntr+1];
	float ite_vec[ntr+1];
	float firstt_vec[ntr+1];
	float ang_vec[ntr+1];
	float angtaper_vec[ntr+1];
	float ang_parcial_vec[ntr+1];
	
	__m256 zero_vec, two_vec, half_vec,one_vec, four_vec, dx_vec, mp_vec, h_vec, v_vec, t0_proc_vec, iip_vec, ip_vec, x_vec, t0_quadrado_vec, ts_sse, tr_sse, t_sse, fnyq_vec;
	__m256 aux_vec0, aux_vec1, aux_vec2, aux_vec3, aux_vec4, aux_vec5, aux_vec6, aux_vec7,aux_vec8,aux_vec9,aux_vec10,aux_vec11,aux_vec12,aux_vec13,aux_vec14,aux_vec15,aux_vec16;
	__m256 p_sse, pmin_sse, nc_vec, dt_vec, three_vec,eigth_vec, itb_sse, nt_vec, fplo_sse, comp_sse, fplo_tmp_sse, ang_sse, twenty_vec, PI_vec, angmax_vec;
	__m256 angtaper_sse, angtaper_tmp_sse, const_180_vec;
	__m256 avx_s;

	/*int calc_indice_amostras = ceil(nt / 4) * 4;*/
	int calc_indice_amostras = ceil(nt / 8) * 8;

	/* Allocate space */
	lowpass 	= alloc2float(nt, nc + 1);
	rtin 		= ealloc1float(nfft);
	rtout 		= ealloc1float(nfft);
	ct 			= ealloc1complex(nf);
	filter 		= alloc1float(nf);
	
	/*
	t_vec= alloc1float(ntr);
	tr_vec= alloc1float(ntr);
	ts_vec= alloc1float(ntr);
	*/
	
	t_vec 		= alloc1float(ntr + 1);
	tr_vec 		= alloc1float(ntr + 1);
	ts_vec 		= alloc1float(ntr + 1);
	geoms_vec 	= alloc1float(ntr + 1);
	obliq_vec 	= alloc1float(ntr + 1);
	pmin_vec 	= alloc1float(ntr + 1);
	//p_vec = alloc1float(ntr+1);
	
	/* Zero all arrays */
	memset((void *) mig[0], 	0, nt*ntr*FSIZE);
	memset((void *) rtin, 		0, nfft * FSIZE);
	//memset((void *) rtout, 		0, nfft * FSIZE);
	memset((void *) filter, 	0, nf * FSIZE);
	memset((void *) lowpass[0], 0, nt * (nc + 1) * FSIZE);
	memset((void *) t_vec, 		0, (ntr + 1) * FSIZE);
	memset((void *) tr_vec, 	0, (ntr + 1) * FSIZE);
	memset((void *) ts_vec, 	0, (ntr + 1) * FSIZE);
	memset((void *) geoms_vec, 	0, (ntr + 1) * FSIZE);
	memset((void *) obliq_vec, 	0, (ntr + 1) * FSIZE);
	memset((void *) pmin_vec, 	0, (ntr + 1) * FSIZE);
	memset((void *) p_vec, 		0, (ntr + 1) * FSIZE);
	memset((void *) fplo_vec, 	0, (ntr + 1) * FSIZE);
	
	/* Preparando variaveis SIMD com constantes */
	one_vec 		= _mm256_set1_ps(1);
	three_vec 		= _mm256_set1_ps(3);
	two_vec 		= _mm256_set1_ps(2);
	four_vec 		= _mm256_set1_ps(4);
	eigth_vec 		= _mm256_set1_ps(8);
	twenty_vec 		= _mm256_set1_ps(20);
	const_180_vec 	= _mm256_set1_ps(180);
	PI_vec 			= _mm256_set1_ps(PI);
	half_vec 		= _mm256_set1_ps(0.5);
	zero_vec 		= _mm256_set1_ps(-0.0);
	nc_vec 			= _mm256_set1_ps(nc);
	dx_vec 			= _mm256_set1_ps(dx);
	dt_vec 			= _mm256_set1_ps(dt);
	nt_vec 			= _mm256_set1_ps(nt);
	angmax_vec 		= _mm256_set1_ps(angmax);
	h_vec 			= _mm256_set1_ps(h);
	/*v_vec 			= _mm256_set1_ps(v);*/
	fnyq_vec 		= _mm256_set1_ps(fnyq);
	
	/* Start the migration process */
	/* Loop over input mid-points first */

	for (imp = 0; imp < ntr; ++imp)	
	{
struct timeval tempo;
		// double tempo_traco;

		tempo = inicia_tempo();

		mp = imp * dx; 
		
		/* Calculate low-pass filtered versions  */
		/* of the data to be used for antialiasing */
		
		//warn("traço de entrada %d/%d\n",imp,ntr);
		
		/*for (it = 0; it < calc_indice_amostras; it += 4)*/
		for (it = 0; it < calc_indice_amostras; it += 8)
		{
			avx_s = _mm256_loadu_ps(&data[imp][it]);
			/*_mm256_store_ps(&rtin[it], avx_s);*/
			_mm256_storeu_ps(&rtin[it], avx_s);
		}

		for (it = calc_indice_amostras; it < nt; ++it){
			rtin[it] = data[imp][it];
		}
				
		/*
		for(it=0; it<nt; ++it){
			rtin[it]=data[imp][it];
		}
		*/ 
			
		for (ifc = 1; ifc < nc + 1; ++ifc)
		{
			memset((void *) rtout, 	0, nfft * FSIZE);
			memset((void *) ct, 	0, nf * FSIZE);
			lpfilt_vec(nfc, nfft, dt, fc[ifc], filter);
			//lpfilt(nfc,nfft,dt,fc[ifc],filter);
			pfarc(1, nfft, rtin, ct);
			for (it = 0; it < nf; ++it) 
			{
				ct[it] = crmul(ct[it], filter[it]);
			}
			pfacr(-1, nfft, ct, rtout);
			/*for (it = 0; it < calc_indice_amostras; it += 8)*/
			for (it = 0; it < calc_indice_amostras; it += 8)
			{
				/*avx_s = _mm256_load_ps(&rtout[it]);*/
				avx_s = _mm256_loadu_ps(&rtout[it]);
				_mm256_storeu_ps(&lowpass[ifc][it], avx_s);
			}

			for (it = calc_indice_amostras; it < nt; ++it)
			{ 
				lowpass[ifc][it] = rtout[it]; 
			}
			/*
			for(it=0; it<nt; ++it){ 
				lowpass[ifc][it]= rtout[it]; 
			}
			*/ 
		}
		
		/* Loop over vertical traveltimes */
		for (it = 0; it < nt; ++it)
		{
			int lx, ux;
			int calculo_indice_abertura;
			//float r,w,g,y,z,q;
			/*float v0,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16;*/
			float ts1,tr1,t0_proc,t0_quadrado;
	
			//tempo = inicia_tempo();
			t0 		= it * dt;
			v 		= vel[imp * dxcdp - firstcdp + cdp_trace_first][it];
			xmax 	= tan((angmax + 10.0) * PI/180.0) * v * t0;
			lx 		= MAX(0,imp - ceil(xmax/dx)); 
			ux 		= MIN(ntr,imp + ceil(xmax/dx));
			//tempo_total = tempo_total + para_tempo(tempo);
			
			t0_proc 					= pow(t0 / 2, 2);
			t0_quadrado 				= t0 * t0;
			/*calculo_indice_abertura 	= ceil(imp / 4) * 4;*/
			/*calculo_indice_abertura 	= ceil(imp / 8) * 8;*/
			calculo_indice_abertura 	= ceil((imp - (lx + 1)) / 8) * 8;
			t0_proc_vec 				= _mm256_set1_ps(t0_proc);
			t0_quadrado_vec 			= _mm256_set1_ps(t0_quadrado);
			mp_vec 						= _mm256_set1_ps(mp);
			v_vec 						= _mm256_set1_ps(v);

			/* for (iip = lx + 1; iip <= calculo_indice_abertura; iip += 4) */
			for (iip = lx + 1; iip <= calculo_indice_abertura; iip += 8)
			{
				/*Calculo do tempo de transito da origem até a subsuperficie*/
				/*iip_vec 		= _mm_set_ps(iip+3, iip+2, iip+1, iip);*/
				iip_vec 		= _mm256_set_ps(iip+7, iip+6, iip+5, iip+4, iip+3, iip+2, iip+1, iip);
				ip_vec 			= _mm256_mul_ps(dx_vec,iip_vec);
				x_vec 			= _mm256_sub_ps(ip_vec,mp_vec);
				aux_vec1 		= _mm256_add_ps(x_vec,h_vec);
				aux_vec2 		= _mm256_div_ps(aux_vec1,v_vec);
				aux_vec3 		= _mm256_mul_ps(aux_vec2,aux_vec2);
				aux_vec4 		= _mm256_add_ps(t0_proc_vec,aux_vec3);
				ts_sse 			= _mm256_sqrt_ps(aux_vec4);
				_mm256_storeu_ps(&ts_vec[iip],ts_sse);
				
				/*Calculo do tempo de transito da subsuperfície até o receptor*/
				aux_vec1 		= _mm256_sub_ps(h_vec,x_vec);
				aux_vec2 		= _mm256_div_ps(aux_vec1,v_vec);
				aux_vec3 		= _mm256_mul_ps(aux_vec2,aux_vec2);
				aux_vec4 		= _mm256_add_ps(t0_proc_vec,aux_vec3);
				tr_sse 			= _mm256_sqrt_ps(aux_vec4);
				_mm256_storeu_ps(&tr_vec[iip],tr_sse);
				
				/*Calculo do tempo de trânsito duplo*/
				t_sse 			= _mm256_add_ps(ts_sse, tr_sse);
				_mm256_storeu_ps(&t_vec[iip],t_sse);
				
				/*Calculo do fator de obliquidade (1)*/
				aux_vec1 		= _mm256_mul_ps(t_sse,v_vec);
				aux_vec2 		= _mm256_div_ps(one_vec,aux_vec1);
				aux_vec3 		= _mm256_sqrt_ps(aux_vec2);
				_mm256_storeu_ps(&geoms_vec[iip],aux_vec3);
				
				/*Calculo do fator de obliquidade (2)*/
				aux_vec0 		= t0_quadrado_vec;
				aux_vec1 		= _mm256_div_ps(aux_vec0,four_vec);
				aux_vec2 		= _mm256_mul_ps(tr_sse,tr_sse);
				aux_vec3 		= _mm256_mul_ps(ts_sse,ts_sse);
				aux_vec4 		= _mm256_mul_ps(tr_sse,ts_sse);
				aux_vec5 		= _mm256_sub_ps(aux_vec2,aux_vec1);
				aux_vec6 		= _mm256_sub_ps(aux_vec3,aux_vec1);
				aux_vec7 		= _mm256_sqrt_ps(aux_vec5);
				aux_vec8 		= _mm256_sqrt_ps(aux_vec6);
				aux_vec9 		= _mm256_div_ps(one_vec,aux_vec4);
				aux_vec10 		= _mm256_div_ps(aux_vec0,_mm256_mul_ps(four_vec,aux_vec4));
				aux_vec11 		= _mm256_mul_ps(aux_vec7,aux_vec8);
				aux_vec12 		= _mm256_mul_ps(aux_vec11,aux_vec9);
				aux_vec13 		= _mm256_sub_ps(aux_vec10,aux_vec12);
				aux_vec14 		= _mm256_add_ps(one_vec,aux_vec13);
				aux_vec15 		= _mm256_mul_ps(half_vec,aux_vec14);
				aux_vec16 		= _mm256_sqrt_ps(aux_vec15);
				_mm256_storeu_ps(&obliq_vec[iip],aux_vec16);
				
				/*Calculo do operador de retardo de migração*/
				
				aux_vec1 		= _mm256_mul_ps(dx_vec,fnyq_vec);
				aux_vec2 		= _mm256_mul_ps(aux_vec1,two_vec);
				pmin_sse 		= _mm256_div_ps(one_vec,aux_vec2);
				_mm256_storeu_ps(&pmin_vec[iip],pmin_sse);
				
				aux_vec0 		= _mm256_add_ps(x_vec,h_vec);
				aux_vec1 		= _mm256_mul_ps(v_vec,v_vec);
				aux_vec2 		= _mm256_mul_ps(aux_vec1,ts_sse);
				aux_vec3 		= _mm256_div_ps(aux_vec0,aux_vec2);
								
				aux_vec4 		= _mm256_sub_ps(x_vec,h_vec);
				aux_vec5 		= _mm256_mul_ps(aux_vec1,tr_sse);
				aux_vec6 		= _mm256_div_ps(aux_vec4,aux_vec5);
				aux_vec7 		= _mm256_add_ps(aux_vec3,aux_vec6);
				p_sse 			= _mm256_andnot_ps(zero_vec,aux_vec7);
				_mm256_storeu_ps(&p_vec[iip],p_sse);
				
				/*Calculo do fplo e fphi*/
				aux_vec1 		= _mm256_mul_ps(nc_vec,pmin_sse);
				aux_vec2 		= _mm256_div_ps(aux_vec1,p_sse);
				fplo_tmp_sse 	= _mm256_floor_ps(aux_vec2);
				//_mm_storeu_ps(&fplo_vec[iip],aux_vec3);
				
				/*comp_sse 		= _mm_cmpgt_ps(p_sse,zero_vec);*/
				comp_sse 		= _mm256_cmp_ps(p_sse, zero_vec, _CMP_GT_OS);
				aux_vec1 		= _mm256_and_ps(fplo_tmp_sse,comp_sse);
				aux_vec2 		= _mm256_andnot_ps(comp_sse,nc_vec);
				fplo_sse 		= _mm256_or_ps(aux_vec1,aux_vec2);
				_mm256_storeu_ps(&fplo_vec[iip],fplo_sse); 

		//		_mm_storeu_ps(&fphi_vec[iip],_mm_add_ps(one_vec,fplo_sse)); 
	
				/*Calculo do intervalo de amostras para interpolar */
				
				aux_vec1 		= _mm256_div_ps(t_sse,dt_vec);
				aux_vec2 		= _mm256_ceil_ps(aux_vec1);
				aux_vec3 		= _mm256_sub_ps(aux_vec2, three_vec);
				itb_sse 		= _mm256_max_ps(aux_vec3, zero_vec);
				_mm256_storeu_ps(&itb_vec[iip],itb_sse);

				aux_vec1 		= _mm256_add_ps(itb_sse,eigth_vec);
				aux_vec2 		= _mm256_min_ps(aux_vec1,nt_vec);
				_mm256_storeu_ps(&ite_vec[iip],aux_vec2);
				
				aux_vec1 		= _mm256_sub_ps(itb_sse,one_vec);
				aux_vec2 		= _mm256_mul_ps(aux_vec1,dt_vec);
				_mm256_storeu_ps(&firstt_vec[iip],aux_vec2);

				/*Calculo do argumento para cômputo da variável ref que é usada para definir os pesos dos filtros*/
				aux_vec1 		= _mm256_mul_ps(nc_vec,pmin_sse);
				_mm256_storeu_ps(&ref_arg_vec[iip],aux_vec1);
			} /* for (iip = lx + 1; iip <= calculo_indice_abertura; iip += 8) */

			for (iip = calculo_indice_abertura + 1; iip <= imp; ++iip)
			{
				ip = iip * dx; 
				x = ip - mp; 
				ts1 = sqrt( pow(t0 / 2, 2) + pow((x + h) / v, 2) );
				tr1 = sqrt( pow(t0 / 2, 2) + pow((h - x) / v, 2) );
				t_vec[iip] = ts1 + tr1;
				tr_vec[iip] = tr1;
				ts_vec[iip] = ts1;
				geoms_vec[iip] = sqrt(1 / (t_vec[iip] * v));
				obliq_vec[iip] = sqrt(.5 * (1 + (t0 * t0 / (4 * ts1 * tr1)) - (1 / (ts1 * tr1)) * sqrt(ts1 * ts1 - t0 * t0 / 4) * sqrt(tr1 * tr1 - t0 * t0 / 4)));
				pmin_vec[iip] = 1/(2*dx*fnyq);
				p_vec[iip] = fabs((x + h) / (pow(v, 2) * ts1) + (x-h) / (pow(v, 2) * tr1));
				if (p_vec[iip] > 0)
					fplo_vec[iip] = floor(nc * pmin_vec[iip] / p_vec[iip]);
				else
					fplo_vec[iip] = nc;
	//			fphi_vec[iip] = fplo_vec[iip]+1;
				//warn("%d %f\n",(int)fphi_vec[iip],fphi_vec[iip]);
				itb_vec[iip] = MAX(ceil(t_vec[iip] / dt)-3, 0);
				ite_vec[iip] = MIN(itb_vec[iip] + 8, nt);
				firstt_vec[iip] = (itb_vec[iip] - 1) * dt;
				ref_arg_vec[iip] = nc * pmin_vec[iip];
			}

			for (iip = imp; iip > lx; --iip)
			{
				//ang_vec[iip]=180.0*fabs(acos(t0/t_vec[iip]))/PI;   
				ang_parcial_vec[iip] = acos(t0 / t_vec[iip]);   
			}
			
			/*for (iip = lx + 1; iip <= calculo_indice_abertura; iip += 4)*/
			for (iip = lx + 1; iip <= calculo_indice_abertura; iip += 8)
			{	
				aux_vec1 			= _mm256_loadu_ps(&ang_parcial_vec[iip]);
				aux_vec2 			= _mm256_andnot_ps(zero_vec, aux_vec1);
				aux_vec3 			= _mm256_mul_ps(aux_vec2, const_180_vec);
				ang_sse 			= _mm256_div_ps(aux_vec3, PI_vec);
				_mm256_storeu_ps(&ang_vec[iip], ang_sse);
				
				//ang_sse = _mm_loadu_ps(&ang_vec[iip]);
				
				aux_vec1 			= _mm256_sub_ps(ang_sse, angmax_vec);
				aux_vec2 			= _mm256_mul_ps(aux_vec1, PI_vec);
				aux_vec3 			= _mm256_div_ps(aux_vec2, twenty_vec);
				angtaper_tmp_sse 	= cos_ps(aux_vec3);
				
				/*comp_sse 			= _mm_cmpgt_ps(ang_sse, angmax_vec);*/
				comp_sse 			= _mm256_cmp_ps(ang_sse, angmax_vec, _CMP_GT_OS);
				aux_vec1 			= _mm256_and_ps(angtaper_tmp_sse,comp_sse);
				aux_vec2 			= _mm256_andnot_ps(comp_sse, one_vec);
				angtaper_sse 		= _mm256_or_ps(aux_vec1, aux_vec2);
				_mm256_storeu_ps(&angtaper_vec[iip], angtaper_sse);
			}
			
			for (iip = calculo_indice_abertura + 1; iip <= imp; ++iip)
			{
				ang_vec[iip] = 180.0 * fabs(acos(t0/t_vec[iip]))/PI;
				if (ang_vec[iip] > angmax)
					angtaper_vec[iip] = cos((ang_vec[iip] - angmax) * PI / 20.0);
				else
					angtaper_vec[iip] = 1.0;
			}
			
			/* loop over output image-points to the left of the midpoint */
			for (iip = imp; iip > lx; --iip)
			{
				/*float ts, tr;*/
				int fplo = 0, fphi = 0/*, fphi2*/;
				float ref, wlo, whi;

				/*ts 	= ts_vec[iip];
				tr 	= tr_vec[iip];*/
				t 	= t_vec[iip];
				if (t >= tmax) break;
				
				geoms = geoms_vec[iip];
				obliq = obliq_vec[iip];
				
//				ang=180.0*fabs(acos(t0/t))/PI;  
//				if(ang<=angmax) angtaper=1.0;
//				if(ang>angmax) angtaper=cos((ang-angmax)*PI/20);

				/*ang = ang_vec[iip];*/
				angtaper = angtaper_vec[iip];

				/* Evaluate migration operator slowness p to determine */
				/* the low-pass filtered trace for antialiasing */
				//pmin=1/(2*dx*fnyq);
				/*pmin = pmin_vec[iip];*/
				p = p_vec[iip];

				fplo = (int)fplo_vec[iip];
				fphi = (int)fplo++;
				
				ref = fmod(ref_arg_vec[iip],p);
				wlo = 1-ref;
				whi = ref;

				itb = itb_vec[iip];
				ite = ite_vec[iip];
				firstt = firstt_vec[iip];
				
				/* Move energy from CMP to CIP */
				if (fplo >= nc)
				{
					if ((ite - itb) == 8)
					{
						avx_s = _mm256_loadu_ps(&lowpass[nc][itb]);
						_mm256_storeu_ps(&datalo[0], avx_s);
						/*
						avx_s = _mm256_loadu_ps(&lowpass[nc][itb + 4]);
						_mm256_store_ps(&datalo[4], avx_s);
						*/
					}
					else
					{					
						for (k = itb; k < ite; ++k){
							datalo[k-itb]=lowpass[nc][k];
						}
					}					
					ints8r(8, dt, firstt, datalo, 0.0, 0.0, 1, &t, &amplo);
					mig[iip][it] += geoms * obliq * angtaper * amplo;
				} 
				else if (fplo < nc)
				{
					if ((ite-itb) == 8)
					{
						avx_s = _mm256_loadu_ps(&lowpass[fplo][itb]);
						_mm256_storeu_ps(&datalo[0],avx_s);
						/*
						avx_s = _mm256_loadu_ps(&lowpass[fplo][itb+4]);
						_mm256_store_ps(&datalo[4],avx_s);
						*/
						avx_s = _mm256_loadu_ps(&lowpass[fphi][itb]);
						_mm256_storeu_ps(&datahi[0],avx_s);
						/*
						avx_s = _mm256_loadu_ps(&lowpass[fphi][itb+4]);
						_mm256_store_ps(&datahi[4],avx_s);
						*/
					}
					else
					{
						for (k = itb; k < ite; ++k)
						{
							datalo[k - itb] = lowpass[fplo][k];
							datahi[k - itb] = lowpass[fphi][k];
						}
					}
					ints8r(8, dt, firstt, datalo, 0.0, 0.0, 1, &t, &amplo);
					ints8r(8, dt, firstt, datahi, 0.0, 0.0, 1, &t, &amphi);
					mig[iip][it] += geoms * obliq * angtaper * (wlo * amplo + whi * amphi);
				}
			}

			/*calculo_indice_abertura = ceil((ux - (imp + 1)) / 4) * 4;*/
			calculo_indice_abertura = ceil((ux - (imp + 1)) / 8) * 8;

			/*for (iip = imp + 1; iip < (imp + 1 + calculo_indice_abertura); iip += 4)*/
			for (iip = imp + 1; iip < (imp + 1 + calculo_indice_abertura); iip += 8)
			{
				/*iip_vec 		= _mm_set_ps(iip + 3, iip + 2, iip + 1, iip);*/
				iip_vec 		= _mm256_set_ps(iip + 7, iip + 6, iip + 5, iip + 4, iip + 3, iip + 2, iip + 1, iip);
				ip_vec 			= _mm256_mul_ps(dx_vec, iip_vec);
				x_vec 			= _mm256_sub_ps(ip_vec, mp_vec);
				aux_vec1 		= _mm256_add_ps(x_vec, h_vec);
				aux_vec2 		= _mm256_div_ps(aux_vec1, v_vec);
				aux_vec3 		= _mm256_mul_ps(aux_vec2, aux_vec2);
				aux_vec4 		= _mm256_add_ps(t0_proc_vec, aux_vec3);
				ts_sse 			= _mm256_sqrt_ps(aux_vec4);
				_mm256_storeu_ps(&ts_vec[iip],ts_sse);
				
				aux_vec1 		= _mm256_sub_ps(h_vec,x_vec);
				aux_vec2 		= _mm256_div_ps(aux_vec1,v_vec);
				aux_vec3 		= _mm256_mul_ps(aux_vec2,aux_vec2);
				aux_vec4 		= _mm256_add_ps(t0_proc_vec,aux_vec3);
				tr_sse 			= _mm256_sqrt_ps(aux_vec4);
				_mm256_storeu_ps(&tr_vec[iip],tr_sse);
				
				t_sse = _mm256_add_ps(ts_sse, tr_sse);
				_mm256_storeu_ps(&t_vec[iip],t_sse);
				
				/*Calculo do fator de obliquidade (1)*/
				aux_vec1 		= _mm256_mul_ps(t_sse,v_vec);
				aux_vec2 		= _mm256_div_ps(one_vec,aux_vec1);
				aux_vec3 		= _mm256_sqrt_ps(aux_vec2);
				_mm256_storeu_ps(&geoms_vec[iip],aux_vec3);
				
				/*Calculo do fator de obliquidade (2)*/
				aux_vec0 		= t0_quadrado_vec;
				aux_vec1 		= _mm256_div_ps(aux_vec0,four_vec);
				aux_vec2 		= _mm256_mul_ps(tr_sse,tr_sse);
				aux_vec3 		= _mm256_mul_ps(ts_sse,ts_sse);
				aux_vec4 		= _mm256_mul_ps(tr_sse,ts_sse);
				
				aux_vec5 		= _mm256_sub_ps(aux_vec2,aux_vec1);
				aux_vec6 		= _mm256_sub_ps(aux_vec3,aux_vec1);
				aux_vec7 		= _mm256_sqrt_ps(aux_vec5);
				aux_vec8 		= _mm256_sqrt_ps(aux_vec6);
				aux_vec9 		= _mm256_div_ps(one_vec,aux_vec4);
				aux_vec10 		= _mm256_div_ps(aux_vec0, _mm256_mul_ps(four_vec, aux_vec4));
				aux_vec11 		= _mm256_mul_ps(aux_vec7,aux_vec8);
				aux_vec12 		= _mm256_mul_ps(aux_vec11,aux_vec9);
				aux_vec13 		= _mm256_sub_ps(aux_vec10,aux_vec12);
				aux_vec14 		= _mm256_add_ps(one_vec,aux_vec13);
				aux_vec15 		= _mm256_mul_ps(half_vec,aux_vec14);
				aux_vec16 		= _mm256_sqrt_ps(aux_vec15);				
				_mm256_storeu_ps(&obliq_vec[iip],aux_vec16);

				/*Calculo do operador de retardo de migração*/
				
				aux_vec1 		= _mm256_mul_ps(dx_vec,fnyq_vec);
				aux_vec2 		= _mm256_mul_ps(aux_vec1,two_vec);
				pmin_sse 		= _mm256_div_ps(one_vec,aux_vec2);
				_mm256_storeu_ps(&pmin_vec[iip],pmin_sse);
				
				aux_vec0 		= _mm256_add_ps(x_vec,h_vec);
				aux_vec1 		= _mm256_mul_ps(v_vec,v_vec);
				aux_vec2 		= _mm256_mul_ps(aux_vec1,ts_sse);
				aux_vec3 		= _mm256_div_ps(aux_vec0,aux_vec2);
				
				aux_vec4 		= _mm256_sub_ps(x_vec,h_vec);
				aux_vec5 		= _mm256_mul_ps(aux_vec1,tr_sse);
				aux_vec6 		= _mm256_div_ps(aux_vec4,aux_vec5);
				aux_vec7 		= _mm256_add_ps(aux_vec3,aux_vec6);
				p_sse 			= _mm256_andnot_ps(zero_vec,aux_vec7);
				_mm256_storeu_ps(&p_vec[iip], p_sse);

				/*Calculo do fplo*/
				aux_vec1 		= _mm256_mul_ps(nc_vec,pmin_sse);
				aux_vec2 		= _mm256_div_ps(aux_vec1,p_sse);
				fplo_tmp_sse 	= _mm256_floor_ps(aux_vec2);
				_mm256_storeu_ps(&fplo_vec[iip],fplo_tmp_sse);

				/*comp_sse 		= _mm_cmpgt_ps(p_sse,zero_vec);*/
				comp_sse 		= _mm256_cmp_ps(p_sse, zero_vec, _CMP_GT_OS);
				aux_vec1 		= _mm256_and_ps(fplo_tmp_sse,comp_sse);
				aux_vec2 		= _mm256_andnot_ps(comp_sse,nc_vec);
				fplo_sse 		= _mm256_or_ps(aux_vec1,aux_vec2);
				_mm256_storeu_ps(&fplo_vec[iip],fplo_sse); 

				/*Calculo do intervalo de amostras para interpolar */
				
				aux_vec1 		= _mm256_div_ps(t_sse,dt_vec);
				aux_vec2 		= _mm256_ceil_ps(aux_vec1);
				aux_vec3 		= _mm256_sub_ps(aux_vec2,three_vec);
				itb_sse 		= _mm256_max_ps(aux_vec3,zero_vec);
				_mm256_storeu_ps(&itb_vec[iip],itb_sse);

				aux_vec1 		= _mm256_add_ps(itb_sse,eigth_vec);
				aux_vec2 		= _mm256_min_ps(aux_vec1,nt_vec);
				_mm256_storeu_ps(&ite_vec[iip],aux_vec2);
				
				aux_vec1 		= _mm256_sub_ps(itb_sse,one_vec);
				aux_vec2 		= _mm256_mul_ps(aux_vec1,dt_vec);
				_mm256_storeu_ps(&firstt_vec[iip],aux_vec2);
				
				/*Calculo do argumento para cômputo da variável ref que é usada para definir os pesos dos filtros*/
				aux_vec1 		= _mm256_mul_ps(nc_vec,pmin_sse);
				_mm256_storeu_ps(&ref_arg_vec[iip],aux_vec1);
			}

			for (iip = (imp + 1 + calculo_indice_abertura); iip < ux; ++iip)
			{
				ip = iip * dx; 
				x = ip - mp; 
				ts1 = sqrt( pow(t0/2,2) + pow((x+h)/v,2) );
				tr1 = sqrt( pow(t0/2,2) + pow((h-x)/v,2) );
				t_vec[iip] = ts1 + tr1;
				tr_vec[iip] = tr1;
				ts_vec[iip] = ts1;
				geoms_vec[iip] = sqrt(1/(t_vec[iip]*v));
				obliq_vec[iip] = sqrt(.5*(1 + (t0*t0/(4*ts1*tr1)) - (1/(ts1*tr1))*sqrt(ts1*ts1 - t0*t0/4)*sqrt(tr1*tr1 - t0*t0/4)));
				pmin_vec[iip] = 1/(2*dx*fnyq);
				p_vec[iip] = fabs((x+h)/(pow(v,2)*ts1) + (x-h)/(pow(v,2)*tr1));
				
				if (p_vec[iip] > 0)
					fplo_vec[iip] = floor(nc * pmin_vec[iip] / p_vec[iip]);
				else
					fplo_vec[iip] = nc;

				itb_vec[iip] = MAX(ceil(t_vec[iip] / dt) - 3, 0);
				ite_vec[iip] = MIN(itb_vec[iip] + 8, nt);
				firstt_vec[iip] = (itb_vec[iip] - 1) * dt;
				ref_arg_vec[iip] = nc * pmin_vec[iip];
			}
			
			for (iip = imp + 1; iip < ux; ++iip)
			{
				ang_parcial_vec[iip] = acos(t0 / t_vec[iip]);   
	
//				ang_vec[iip]=180.0*fabs(acos(t0/t_vec[iip]))/PI;   
//				v1 = (ang_vec[iip]-angmax);
//				v2 = v1*PI;
//				v3 = v2/20.0;
//				v4 = cos(v3);
				
//				angtaper_vec[iip]=v4;

			}
			
			/*for (iip = imp + 1; iip < (imp + 1 + calculo_indice_abertura); iip += 4)*/
			for (iip = imp + 1; iip < (imp + 1 + calculo_indice_abertura); iip += 8)
			{
				aux_vec1 		= _mm256_loadu_ps(&ang_parcial_vec[iip]);
				aux_vec2 		= _mm256_andnot_ps(zero_vec,aux_vec1);
				aux_vec3 		= _mm256_mul_ps(aux_vec2,const_180_vec);
				ang_sse 		= _mm256_div_ps(aux_vec3,PI_vec);
				_mm256_storeu_ps(&ang_vec[iip],ang_sse);
								
				//ang_sse = _mm_loadu_ps(&ang_vec[iip]);
				aux_vec1 		= _mm256_sub_ps(ang_sse,angmax_vec);
				aux_vec2 		= _mm256_mul_ps(aux_vec1,PI_vec);
				aux_vec3 		= _mm256_div_ps(aux_vec2,twenty_vec);
				angtaper_tmp_sse = cos_ps(aux_vec3);
				
				/*comp_sse 		= _mm_cmpgt_ps(ang_sse,angmax_vec);*/
				comp_sse 		= _mm256_cmp_ps(ang_sse, angmax_vec, _CMP_GT_OS);
				aux_vec1 		= _mm256_and_ps(angtaper_tmp_sse, comp_sse);
				aux_vec2 		= _mm256_andnot_ps(comp_sse, one_vec);
				angtaper_sse 	= _mm256_or_ps(aux_vec1, aux_vec2);
				_mm256_storeu_ps(&angtaper_vec[iip], angtaper_sse);
			}
			
			for (iip = (imp + 1 + calculo_indice_abertura); iip < ux; ++iip)
			{
				ang_vec[iip] = 180.0 * fabs(acos(t0 / t_vec[iip])) / PI;   
				
				if (ang_vec[iip] > angmax) 
					angtaper_vec[iip] = cos((ang_vec[iip] - angmax) * PI / 20.0);
				else
					angtaper_vec[iip] = 1.0;
			}

			/* loop over output image-points to the right of the midpoint */
			for (iip = imp + 1; iip < ux; ++iip)
			{
				/*float ts,tr;*/
				int fplo = 0, fphi;
				float ref, wlo, whi;
				
				/*ts = ts_vec[iip];*/
				/*tr = tr_vec[iip];*/
				t = t_vec[iip];
				
				if (t >= tmax) break;
				geoms = geoms_vec[iip];
				obliq = obliq_vec[iip];

//				tempo = inicia_tempo();
				//ang=180.0*fabs(acos(t0/t))/PI;   
				/*ang = ang_vec[iip];*/
				angtaper = angtaper_vec[iip];

				//tempo = inicia_tempo();
				//tempo_total = tempo_total + para_tempo(tempo);

				/* Evaluate migration operator slowness p to determine the  */
				/* low-pass filtered trace for antialiasing */
				/*pmin = pmin_vec[iip];*/
				p = p_vec[iip];
				fplo = fplo_vec[iip];
				fphi = fplo + 1;
	
				ref = fmod(ref_arg_vec[iip], p);
				
				wlo = 1 - ref;
				whi = ref;
				
				itb = itb_vec[iip];
				ite = ite_vec[iip];
				firstt = firstt_vec[iip];

				/* Move energy from CMP to CIP */
				if (fplo >= nc)
				{
					if ((ite-itb) == 8)
					{
						avx_s = _mm256_loadu_ps(&lowpass[nc][itb]);
						_mm256_storeu_ps(&datalo[0], avx_s);
						/*
						avx_s = _mm256_loadu_ps(&lowpass[nc][itb + 4]);
						_mm256_store_ps(&datalo[4], avx_s);
						*/
					}
					else
					{
						for (k = itb; k < ite; ++k)
						{
							datalo[k-itb] = lowpass[nc][k];
						}
					}
					ints8r(8, dt, firstt, datalo, 0.0, 0.0, 1, &t, &amplo);
					mig[iip][it] += geoms * obliq * angtaper * amplo;
				} 
				else if (fplo < nc) 
				{
					if ((ite-itb) == 8)
					{
						avx_s = _mm256_loadu_ps(&lowpass[fplo][itb]);
						_mm256_storeu_ps(&datalo[0], avx_s);
						/*
						avx_s = _mm256_loadu_ps(&lowpass[fplo][itb + 4]);
						_mm256_store_ps(&datalo[4], avx_s);
						*/
						avx_s = _mm256_loadu_ps(&lowpass[fphi][itb]);
						_mm256_storeu_ps(&datahi[0], avx_s);
						/*
						avx_s = _mm256_loadu_ps(&lowpass[fphi][itb + 4]);
						_mm256_store_ps(&datahi[4], avx_s);
						*/
					}
					else
					{
						for (k = itb; k < ite; ++k)
						{
							datalo[k - itb] = lowpass[fplo][k];
							datahi[k - itb] = lowpass[fphi][k];
						}
					}
					ints8r(8, dt, firstt, datalo, 0.0, 0.0, 1, &t, &amplo);
					ints8r(8, dt, firstt, datahi, 0.0, 0.0, 1, &t, &amphi);
					mig[iip][it] += geoms * obliq * angtaper * (wlo * amplo + whi * amphi);
				}
			} /* for (iip = imp + 1; iip < ux; ++iip) */ /* loop over output image-points to the right of the midpoint */
		} /* for (it = 0; it < nt; ++it) */ /* Loop over vertical traveltimes */
tempo_offset += para_tempo(tempo);
	} /* for (imp = 0; imp < ntr; ++imp) */ /* Loop over input mid-points first */
} /* void  migkt2d (...) */

void lpfilt_vec(int nfc, int nfft, float dt, float fhi, float *filter)
{
	int calculo_indice_nf;
	int calculo_indice_nfc;
	int i,j;  /* counters */
	int nf;   /* Number of frequencies (including Nyquist) */
	float onfft;  /* reciprocal of nfft */
	float fn; /* Nyquist frequency */
	float df; /* frequency interval */
	float dw; /* frequency interval in radians */
	float whi;/* cut-frequency in radians */
	float whi_conv;
	//float dobro_nfc;
	/*__m128 whi_conv_vec, dw_vec, indice_nf_vec, indice_nfc_vec, aux, aux2;*/
	__m256 whi_conv_vec, dw_vec, indice_nf_vec, indice_nfc_vec, aux, aux2;
	
	nf = nfft / 2 + 1;
	onfft = 1.0 / nfft;
	fn = 1.0 / (2 * dt);
	df = onfft / dt;
	whi = fhi * PI / fn;
	dw = df * PI / fn;
	whi_conv = whi / PI;
	//dobro_nfc = nfc*2;
	/*calculo_indice_nf = ceil(nf / 4) * 4;*/
	calculo_indice_nf = ceil(nf / 8) * 8;
	/*calculo_indice_nfc = ceil(nfc / 4) * 4;*/
	calculo_indice_nfc = ceil(nfc / 8) * 8;

	/*float c_vec[nfc] __attribute__((aligned(16)));*/
	float c_vec[nfc] __attribute__((aligned(32)));
	//float c_vec2[nfc] __attribute__((aligned(16)));;	
	/*float w_vec[nf] __attribute__((aligned(16)));*/
	float w_vec[nf] __attribute__((aligned(32)));
	/*float parcial[nfc] __attribute__((aligned(16)));*/
	float parcial[nfc] __attribute__((aligned(32)));	
	/*float indice_nf[nf] __attribute__((aligned(16)));*/
	float indice_nf[nf] __attribute__((aligned(32)));
	/*float indice_nfc[nfc] __attribute__((aligned(16)));*/
	float indice_nfc[nfc] __attribute__((aligned(32)));
	//float v1, v2;
	//float v3, v4;
	
	for (i = 0; i < nf; i++)
		indice_nf[i] = i;
	for (i = 0; i < nfc; i++)
		indice_nfc[i] = i;
	
	/*
	for(j=1; j<calculo_indice_nfc; j+=4)
	{
		//c_vec[j]= sin(whi*j) * sin(PI*j/nfc) * (2*nfc) / (PI*PI*j*j);
		aux1 = _sin_ps(_mm_mul_ps(whi_vec,indice_nfc_vec));
		aux2 = _sin_ps(_mm_div_ps(_mm_mul_ps(PI_vec,indice_nfc_vec),nfc_vec));
		aux3 = _mm_mul_ps(_pow_ps(PI_vec),_pow_ps(indice_nfc_vec));
		resultado = _mm_mul_ps(aux1,_mm_mul_ps(aux2,_mm_div_ps(dobro_nfc_vec,aux3)));
		
		_mm_store_ps(&c_vec[j], aux);
	}
	
	for(j=calculo_indice_nfc; j<nfc; j++)
	{
		c_vec[j]= sin(whi*j)*sin(PI*j/nfc)*2*nfc/(PI*PI*j*j);
	}
	*/
	
	//tempo = inicia_tempo();
	for (j = 1; j < nfc; ++j)
	{
		c_vec[j] = sin(whi * j) * sin(PI * j / nfc) * 2 * nfc / (PI * PI * j * j);
	}
	//tempo_total = tempo_total + para_tempo(tempo);
	
	whi_conv_vec = _mm256_set1_ps(whi_conv);
	dw_vec = _mm256_set1_ps(dw);

	/*for (i = 0; i < calculo_indice_nf; i += 4)*/
	for (i = 0; i < calculo_indice_nf; i += 8)
	{
		_mm256_storeu_ps(&filter[i], whi_conv_vec);
		indice_nf_vec = _mm256_loadu_ps(&indice_nf[i]);
		_mm256_storeu_ps(&w_vec[i], _mm256_mul_ps(indice_nf_vec,dw_vec));
		
	}
	
	for (i = calculo_indice_nf; i < nf; ++i)
	{
		filter[i] = whi_conv;
		w_vec[i] = i * dw;
	}

	for (i = 0; i < nf; ++i)
	{		
		aux2 = _mm256_set1_ps(w_vec[i]);
		/*for (j = 1; j < calculo_indice_nfc; j += 4)*/
		for (j = 1; j < calculo_indice_nfc; j += 8)
		{
			aux = _mm256_loadu_ps(&c_vec[j]);
			indice_nfc_vec = _mm256_loadu_ps(&indice_nfc[j]);
			_mm256_storeu_ps(&parcial[j],_mm256_mul_ps(aux, cos_ps(_mm256_mul_ps(indice_nfc_vec, aux2))));
		}
		
		for (j = calculo_indice_nfc; j < nfc; ++j)
			parcial[j] = c_vec[j] * cos(j * w_vec[i]);
	
		for (j = 1; j < nfc; ++j)
			filter[i] += parcial[j];
	}
} /*void lpfilt_vec(int nfc, int nfft, float dt, float fhi, float *filter)*/

void
lpfilt(int nfc, int nfft, float dt, float fhi, float *filter)
/*******************************************************************************
lpfilt -- low-pass filter using Lanczos Smoothing
	(R.W. Hamming:"Digital Filtering",1977)
****************************************************************************
Input:
nfc	number of Fourier coefficients to approximate ideal filter
nfft	number of points in the fft
dt	time sampling interval
fhi	cut-frequency

Output:
filter  array[nf] of filter values
*****************************************************************************
Notes: Filter is to be applied in the frequency domain
*****************************************************************************
Author: CWP: Carlos Pacheco   2006
*****************************************************************************/
{
	int i,j;  /* counters */
	int nf;   /* Number of frequencies (including Nyquist) */
	float onfft;  /* reciprocal of nfft */
	float fn; /* Nyquist frequency */
	float df; /* frequency interval */
	float dw; /* frequency interval in radians */
	float whi;/* cut-frequency in radians */
	float w;  /* radian frequency */

	nf= nfft/2 + 1;
	onfft=1.0/nfft;
	fn=1.0/(2*dt);
	df=onfft/dt;
	whi=fhi*PI/fn;
	dw=df*PI/fn;

	for(i=0; i<nf; ++i){
		filter[i]= whi/PI;
		w=i*dw;

		for(j=1; j<nfc; ++j){
			float c= sin(whi*j)*sin(PI*j/nfc)*2*nfc/(PI*PI*j*j);
			filter[i] +=c*cos(j*w);
		}
	}
}

struct timeval inicia_tempo()
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv;	
}

double para_tempo(struct timeval tv1)
{
	struct timeval tv2;
	float tempo;
	gettimeofday(&tv2, NULL);
	tempo = (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +  (double) (tv2.tv_sec - tv1.tv_sec);
	return tempo;
}

void gravar_tempo_offset(int oldoffset, double tempo, char *arq)
{
	FILE *arquivo;
	arquivo = fopen(arq,"a");
	fprintf(arquivo,"%f;",tempo);
}

void gravar_dados(float **mig, int ntr, int nt, int oldoffset, char *arq)
{
	FILE *arquivo;
	warn("Offset=%d, ntr=%d, nt=%d", oldoffset, ntr, nt);
	arquivo = fopen(arq,"a");
	int ix, it;
	for (ix=0; ix<ntr; ++ix) 
	{
		fprintf(arquivo,"Inicio Traco %d - Offset %d\n",ix,oldoffset);
		for (it=0; it<nt; ++it) {
			fprintf(arquivo,"%f ",mig[ix][it]);
		}
		fprintf(arquivo,"\nFim Traco %d - Offset %d\n",ix,oldoffset);
		
	}
	fclose(arquivo);	
}
