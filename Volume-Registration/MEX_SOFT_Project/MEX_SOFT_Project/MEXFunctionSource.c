

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "FFTWLibraries/fftw3.h"
//#include "Include/csecond.h"
#include <time.h>
#include "Include/makeweights.h"
#include "Include/so3_correlate_fftw.h"
#include "Include/soft_fftw.h"

#include "Include/s2_cospmls.h"
#include "Include/s2_legendreTransforms.h"
#include "Include/s2_semi_memo.h"

#include "mex.h"                                 /*******   This is MATLAB business   ********/

#define NORM( x ) ( (x[0])*(x[0]) + (x[1])*(x[1]) )

// To enable this MEX file computation change the configuration from  .exe to .dll 

void ComputeCorrelatioOnSphere(double* dSignalRealPart, double* dSignalImagPart, double* dPatternRealPart, double* dPatternImagPart, double bandwidthForward, double bandwidthInverse, double degrees, double* alpha, double* beta, double* gamma)
{
	// Function extracted from the SO(3) source with minimal changes that are required for the correlation to work from external data source such as passed input arguments from MATLAB
	int i;
	int n, bwIn, bwOut, degLim;
	double tstart, tstop;
	fftw_complex *workspace1, *workspace2;
	double *workspace3;
	double *sigR, *sigI;
	double *sigCoefR, *sigCoefI;
	double *patCoefR, *patCoefI;
	fftw_complex *so3Sig, *so3Coef;
	fftw_plan p1;
	int na[2], inembed[2], onembed[2];
	int rank, howmany, istride, idist, ostride, odist;
	int tmp, maxloc, ii, jj, kk;
	double maxval, tmpval;
	double *weights;
	double *seminaive_naive_tablespace;
	double **seminaive_naive_table;
	fftw_plan dctPlan, fftPlan;
	int howmany_rank;
	fftw_iodim dims[1], howmany_dims[1];

	bwIn = bandwidthForward;
	bwOut = bandwidthInverse;
	degLim = degrees;


	n = 2 * bwIn;

	sigR = (double *)calloc(n * n, sizeof(double));
	sigI = (double *)calloc(n * n, sizeof(double));
	so3Sig = fftw_malloc(sizeof(fftw_complex) * (8 * bwOut*bwOut*bwOut));
	workspace1 = fftw_malloc(sizeof(fftw_complex) * (8 * bwOut*bwOut*bwOut));
	workspace2 = fftw_malloc(sizeof(fftw_complex) * ((14 * bwIn*bwIn) + (48 * bwIn)));
	workspace3 = (double *)malloc(sizeof(double) * (12 * n + n*bwIn));
	sigCoefR = (double *)malloc(sizeof(double) * bwIn * bwIn);
	sigCoefI = (double *)malloc(sizeof(double) * bwIn * bwIn);
	patCoefR = (double *)malloc(sizeof(double) * bwIn * bwIn);
	patCoefI = (double *)malloc(sizeof(double) * bwIn * bwIn);
	so3Coef = fftw_malloc(sizeof(fftw_complex) * ((4 * bwOut*bwOut*bwOut - bwOut) / 3));


	seminaive_naive_tablespace =
		(double *)malloc(sizeof(double) *
		(Reduced_Naive_TableSize(bwIn, bwIn) +
		Reduced_SpharmonicTableSize(bwIn, bwIn)));

	weights = (double *)malloc(sizeof(double) * (4 * bwIn));

	/****
	At this point, check to see if all the memory has been
	allocated. If it has not, there's no point in going further.
	****/

	if ((seminaive_naive_tablespace == NULL) || (weights == NULL) ||
		(sigR == NULL) || (sigI == NULL) ||
		(so3Coef == NULL) ||
		(workspace1 == NULL) || (workspace2 == NULL) ||
		(workspace3 == NULL) ||
		(sigCoefR == NULL) || (sigCoefI == NULL) ||
		(patCoefR == NULL) || (patCoefI == NULL) ||
		(so3Sig == NULL))
	{
		perror("Error in allocating memory");
		exit(1);
	}

	/* create fftw plans for the S^2 transforms */
	/* first for the dct */
	dctPlan = fftw_plan_r2r_1d(2 * bwIn, weights, workspace3,
		FFTW_REDFT10, FFTW_ESTIMATE);

	/* now for the fft */
	/*
	IMPORTANT NOTE!!! READ THIS!!!

	Now to make the fft plans.

	Please note that the planning-rigor flag *must be* FFTW_ESTIMATE!
	Why? Well, to try to keep things simple. I am using some of the
	pointers to arrays in rotateFct's arguments in the fftw-planning
	routines. If the planning-rigor is *not* FFTW_ESTIMATE, then
	the arrays will be written over during the planning stage.

	Therefore, unless you are really really sure you know what
	you're doing, keep the rigor as FFTW_ESTIMATE !!!
	*/

	/*
	fftw "preamble" ;
	note  that this places in the transposed array
	*/

	rank = 1;
	dims[0].n = 2 * bwIn;
	dims[0].is = 1;
	dims[0].os = 2 * bwIn;
	howmany_rank = 1;
	howmany_dims[0].n = 2 * bwIn;
	howmany_dims[0].is = 2 * bwIn;
	howmany_dims[0].os = 1;

	fftPlan = fftw_plan_guru_split_dft(rank, dims,
		howmany_rank, howmany_dims,
		sigR, sigI,
		(double *)workspace2,
		(double *)workspace2 + (n*n),
		FFTW_ESTIMATE);

	/* create plan for inverse SO(3) transform */
	n = 2 * bwOut;
	howmany = n*n;
	idist = n;
	odist = n;
	rank = 2;
	inembed[0] = n;
	inembed[1] = n*n;
	onembed[0] = n;
	onembed[1] = n*n;
	istride = 1;
	ostride = 1;
	na[0] = 1;
	na[1] = n;

	p1 = fftw_plan_many_dft(rank, na, howmany,
		workspace1, inembed,
		istride, idist,
		so3Sig, onembed,
		ostride, odist,
		FFTW_FORWARD, FFTW_ESTIMATE);


	fprintf(stdout, "Generating seminaive_naive tables...\n");
	seminaive_naive_table = SemiNaive_Naive_Pml_Table(bwIn, bwIn,
		seminaive_naive_tablespace,
		(double *)workspace2);


	/* make quadrature weights for the S^2 transform */
	makeweights(bwIn, weights);

	n = 2 * bwIn;
	//printf("Reading in signal file\n");
	/* read in SIGNAL samples */
	/* first the signal */
	/*fp = fopen(argv[1], "r");
	for (i = 0; i < n * n; i++)
	{
		fscanf(fp, "%lf", sigR + i);
		fscanf(fp, "%lf", sigI + i);
	}
	fclose(fp);*/

	mexPrintf("now taking spherical transform of signal\n");
	FST_semi_memo(dSignalRealPart, dSignalImagPart,
		sigCoefR, sigCoefI,
		bwIn, seminaive_naive_table,
		(double *)workspace2, 0, bwIn,
		&dctPlan, &fftPlan,
		weights);

	//printf("Reading in pattern file\n");
	///* read in SIGNAL samples */
	///* first the signal */
	//fp = fopen(argv[2], "r");
	//for (i = 0; i < n * n; i++)
	//{
	//	fscanf(fp, "%lf", sigR + i);
	//	fscanf(fp, "%lf", sigI + i);
	//}
	//fclose(fp);

	mexPrintf("now taking spherical transform of pattern\n");
	FST_semi_memo(dPatternRealPart, dPatternImagPart,
		patCoefR, patCoefI,
		bwIn, seminaive_naive_table,
		(double *)workspace2, 0, bwIn,
		&dctPlan, &fftPlan,
		weights);

	mexPrintf("freeing seminaive_naive_table and seminaive_naive_tablespace\n");

	free(seminaive_naive_table);
	free(seminaive_naive_tablespace);


	mexPrintf("about to combine coefficients\n");

	/* combine coefficients */
	//tstart = clock() / CLOCKS_PER_SEC; //  csecond() ;
	so3CombineCoef_fftw(bwIn, bwOut, degLim,
		sigCoefR, sigCoefI,
		patCoefR, patCoefI,
		so3Coef);
	//tstop = clock() / CLOCKS_PER_SEC; //  csecond();
	//fprintf(stderr, "combine time \t = %.4e\n", tstop - tstart);

	mexPrintf("about to inverse so(3) transform\n");

	//tstart = clock() / CLOCKS_PER_SEC; //  csecond();
	/* now inverse so(3) */
	Inverse_SO3_Naive_fftw(bwOut,
		so3Coef,
		so3Sig,
		workspace1,
		workspace2,
		workspace3,
		&p1,
		0);
	//tstop = clock() / CLOCKS_PER_SEC; // csecond();
	mexPrintf("finished inverse so(3) transform\n");
	//fprintf(stderr, "inverse so(3) time \t = %.4e\n", tstop - tstart);


	/* now find max value */
	maxval = 0.0;
	maxloc = 0;
	for (i = 0; i < 8 * bwOut*bwOut*bwOut; i++)
	{
		/*
		if (so3Sig[i][0] >= maxval)
		{
		maxval = so3Sig[i][0];
		maxloc = i ;
		}
		*/
		tmpval = NORM(so3Sig[i]);
		if (tmpval > maxval)
		{
			maxval = tmpval;
			maxloc = i;
		}

	}

	ii = floor(maxloc / (4.*bwOut*bwOut));
	tmp = maxloc - (ii*4.*bwOut*bwOut);
	jj = floor(tmp / (2.*bwOut));
	tmp = maxloc - (ii * 4 * bwOut*bwOut) - jj*(2 * bwOut);
	kk = tmp;

	mexPrintf("ii = %d\tjj = %d\tkk = %d\n", ii, jj, kk);

	alpha[0] = M_PI*jj / ((double)bwOut);
	beta[0] = M_PI*(2 * ii + 1) / (4.*bwOut);
	gamma[0] = M_PI*kk / ((double)bwOut);
	mexPrintf("alpha = %f\nbeta = %f\ngamma = %f\n",
		*alpha,
		*beta,
		*gamma);


	///* now save data -> just the real part because the
	//imaginary parts should all be 0 ... right?!? */
	//if (argc == 7)
	//{
	//	printf("about to save data\n");
	//	fp = fopen(argv[6], "w");
	//	for (i = 0; i < 8 * bwOut*bwOut*bwOut; i++)
	//		fprintf(fp, "%.16f\n", so3Sig[i][0]);
	//	fclose(fp);
	//}


	fftw_destroy_plan(p1);
	fftw_destroy_plan(fftPlan);
	fftw_destroy_plan(dctPlan);
	

	free(weights);
	fftw_free(so3Coef);
	free(patCoefI);
	free(patCoefR);
	free(sigCoefI);
	free(sigCoefR);
	free(workspace3);
	fftw_free(workspace2);
	fftw_free(workspace1);
	fftw_free(so3Sig);
	free(sigI);
	free(sigR);

}

void HandleCorrelationOnSphere(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// Extract signal and pattern information from the MATLAB Mex passed arguments and send for actual Correlation Processing
	double* dSignalRealPart;
	double*	dSignalImagPart;
	double* dPatternRealPart;
	double* dPatternImagPart;    // Input Arrays
	double bandwidthForward, bandwidthInverse, degrees;                                 // Input scalars
	double* alpha = 0;
	double* beta = 0;
	double* gamma = 0;														// Output scalars		

	/* create a pointer to the real data in the input matrix  */
	dSignalRealPart = mxGetPr(prhs[1]);
	dSignalImagPart = mxGetPr(prhs[2]);
	dPatternRealPart = mxGetPr(prhs[3]);
	dPatternImagPart = mxGetPr(prhs[4]);

	/* get the value of the scalar input  */
	bandwidthForward = mxGetScalar(prhs[5]);
	bandwidthInverse = mxGetScalar(prhs[6]);
	degrees = mxGetScalar(prhs[7]);

	/* Assign a pointer to the output */
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);


	alpha = mxGetPr(plhs[0]);
	beta = mxGetPr(plhs[1]);
	gamma = mxGetPr(plhs[2]);

	/* call the computational routine */
	ComputeCorrelatioOnSphere(dSignalRealPart, dSignalImagPart, dPatternRealPart, dPatternImagPart, bandwidthForward, bandwidthInverse, degrees, alpha, beta, gamma);
}

void ComputeRotationOnSphere(double* dSignalRealPart, double* dSignalImagPart, double bandwidthForward, double bandwidthInverse, double degrees, double alpha_given, double beta_given, double gamma_given, double* dPatternRealPart, double* dPatternImagPart)
{
	//FILE *fp;
	int i;
	int bwIn, bwOut, degOut;
	double alpha, beta, gamma;
	double *sigInR, *sigInI, *sigOutR, *sigOutI;
	double *scratch;
	double tstart, tstop;
	double *seminaive_naive_tablespace, *trans_seminaive_naive_tablespace2;
	double *seminaive_naive_tablespace2;
	double **seminaive_naive_table2, **seminaive_naive_table;
	double **trans_seminaive_naive_table2;

	
	bwIn = (int)bandwidthForward;
	bwOut = (int)bandwidthInverse;
	degOut = (int)degrees;
	alpha = alpha_given;
	beta = beta_given;
	gamma = gamma_given;

	sigInR = (double *)malloc(sizeof(double)*(4 * bwIn*bwIn));
	sigInI = (double *)malloc(sizeof(double)*(4 * bwIn*bwIn));
	sigOutR = (double *)malloc(sizeof(double)*(4 * bwOut*bwOut));
	sigOutI = (double *)malloc(sizeof(double)*(4 * bwOut*bwOut));

	if (bwOut > bwIn)
		scratch = (double *)malloc(sizeof(double)*((14 * bwOut*bwOut) + (52 * bwOut) + (2 * bwIn)));
	else
		scratch = (double *)malloc(sizeof(double)*((14 * bwIn*bwIn) + (52 * bwIn) + (2 * bwIn)));

	seminaive_naive_tablespace =
		(double *)malloc(sizeof(double) *
		(Reduced_Naive_TableSize(bwIn, bwIn) +
		Reduced_SpharmonicTableSize(bwIn, bwIn)));

	trans_seminaive_naive_tablespace2 =
		(double *)malloc(sizeof(double) *
		(Reduced_Naive_TableSize(bwOut, bwOut) +
		Reduced_SpharmonicTableSize(bwOut, bwOut)));

	seminaive_naive_tablespace2 =
		(double *)malloc(sizeof(double) *
		(Reduced_Naive_TableSize(bwOut, bwOut) +
		Reduced_SpharmonicTableSize(bwOut, bwOut)));



	/****
	At this point, check to see if all the memory has been
	allocated. If it has not, there's no point in going further.
	****/

	if ((scratch == NULL) ||
		(sigInR == NULL) || (sigInI == NULL) ||
		(sigOutR == NULL) || (sigOutI == NULL) ||
		(seminaive_naive_tablespace == NULL) ||
		(trans_seminaive_naive_tablespace2 == NULL))
	{
		perror("Error in allocating memory");
		exit(1);
	}


	mexPrintf(stdout, "Generating seminaive_naive tables...\n");
	seminaive_naive_table = SemiNaive_Naive_Pml_Table(bwIn, bwIn,
		seminaive_naive_tablespace,
		scratch);


	mexPrintf(stdout, "Generating seminaive_naive tables...\n");
	seminaive_naive_table2 = SemiNaive_Naive_Pml_Table(bwOut, bwOut,
		seminaive_naive_tablespace2,
		scratch);


	mexPrintf(stdout, "Generating trans_seminaive_naive tables...\n");
	trans_seminaive_naive_table2 =
		Transpose_SemiNaive_Naive_Pml_Table(seminaive_naive_table2,
		bwOut, bwOut,
		trans_seminaive_naive_tablespace2,
		scratch);

	//fprintf(stdout, "reading in signal ...\n");

	///* read in signal */
	//fp = fopen(argv[7], "r");
	//for (i = 0; i < (4 * bwIn*bwIn); i++)
	//{
	//	fscanf(fp, "%lf", sigInR + i);
	//	fscanf(fp, "%lf", sigInI + i);
	//}
	//fclose(fp);

	mexPrintf(stdout, "about to rotate ...\n");
	
	//tstart = clock() / CLOCKS_PER_SEC; // clock() / CLOCKS_PER_SEC; // csecond();
	rotateFctFFTW(bwIn, bwOut, degOut,
		dSignalRealPart, dSignalImagPart,
		sigOutR, sigOutI,
		alpha, beta, gamma,
		scratch,
		seminaive_naive_table,
		trans_seminaive_naive_table2);

	//tstop = clock() / CLOCKS_PER_SEC; // csecond();
	mexPrintf(stdout, "finished rotating ...\n");
	//fprintf(stdout, "rotation time \t = %.4e\n", tstop - tstart);

	///* write out rotated signal */
	//fp = fopen(argv[8], "w");
	//for (i = 0; i < (4 * bwOut*bwOut); i++)
	//{
	//	fprintf(fp, "%.15f\n%.15f\n", sigOutR[i], sigOutI[i]);
	//}
	//fclose(fp);

	memcpy(dPatternRealPart, sigOutR, sizeof(double)*(4 * bwOut*bwOut));
	memcpy(dPatternImagPart, sigOutI, sizeof(double)*(4 * bwOut*bwOut));


	mexPrintf(stdout, "finished writing ...\n");

	free(trans_seminaive_naive_table2);
	free(seminaive_naive_table2);
	free(seminaive_naive_table);
	free(seminaive_naive_tablespace2);
	free(trans_seminaive_naive_tablespace2);
	free(seminaive_naive_tablespace);

	free(scratch);
	free(sigOutI);
	free(sigOutR);
	free(sigInI);
	free(sigInR);
}

void HandleRotationOnSphere(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double bandwidthForward, bandwidthInverse, degrees;                                 // Input scalars
	double alpha , beta, gamma;

	double* dSignalRealPart;
	double*	dSignalImagPart;

	double* dPatternRealPart;                      // Output Arrays
	double* dPatternImagPart;						// Output Arrays

	/* get the value of the scalar input  */
	bandwidthForward = mxGetScalar(prhs[1]);
	bandwidthInverse = mxGetScalar(prhs[2]);
	degrees = mxGetScalar(prhs[3]);
	alpha = mxGetScalar(prhs[4]);
	beta = mxGetScalar(prhs[5]);
	gamma = mxGetScalar(prhs[6]);

	/* create a pointer to the real data in the input matrix  */
	dSignalRealPart = mxGetPr(prhs[7]);
	dSignalImagPart = mxGetPr(prhs[8]);

	/* Assign a pointer to the output */
	plhs[0] = mxCreateDoubleMatrix(1, 4 * bandwidthInverse*bandwidthInverse, mxREAL);         //   This should be the output size 
	plhs[1] = mxCreateDoubleMatrix(1, 4 * bandwidthInverse*bandwidthInverse, mxREAL);		//   This should be the output size 
	
	dPatternRealPart = mxGetPr(plhs[0]);
	dPatternImagPart = mxGetPr(plhs[1]);

	ComputeRotationOnSphere(dSignalRealPart, dSignalImagPart, bandwidthForward, bandwidthInverse, degrees, alpha, beta, gamma, dPatternRealPart, dPatternImagPart);

}



/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int IsThisRotation = mxGetScalar(prhs[0]);             // First value determines if we are rotating or estimating, since MEX project allows only one mex file 

	if (0 == IsThisRotation)   
		// This is estimation from two signals on a sphere, sytax to run in C-Source SO(3) test_soft_fftw_correlate2.c is:
		// randomS2sigA_bw8.dat randomS2sig_bw8.dat 8 8 7   ==   myRotated.dat randomS2sig_bw8.dat 8 8 7
	{
		/* check for proper number of arguments */
		if (nrhs != 8) {
			mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "Seven inputs required.");
		}
		if (nlhs != 3) {
			mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs", "Three output required.");
		}
		HandleCorrelationOnSphere( nlhs, plhs, nrhs, prhs);
	}
	else						
		// This is rotation of a single signal on a sphere, sytax in C-Source SO(3) test_s2_rotate_fftw.c is:
		// 8 8 7 0.392699 1.07922 0.785398 randomS2sig_bw8.dat  myRotated.dat 
	{
		/* check for proper number of arguments */
		if (nrhs != 9) {
			mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "Seven inputs required.");
		}
		if (nlhs != 2) {
			mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs", "Two outputs as returned array required.");
		}
		HandleRotationOnSphere(nlhs, plhs, nrhs, prhs);
	}

}