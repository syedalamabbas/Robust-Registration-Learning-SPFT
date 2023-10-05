///*==========================================================
//* Similar to arrayProduct.c - example in MATLAB External Interfaces
//*
//* Multiplies an input scalar (multiplier)
//* times a 1xN matrix (inMatrix)
//* and outputs a 1xN matrix (outMatrix)
//*
//* The calling syntax is:
//*
//*		outMatrix = arrayProduct(multiplier, inMatrix)
//*
//* This is a MEX-file for MATLAB.
//* Copyright 2007-2010 The MathWorks, Inc.
//*
//*========================================================*/
///* $Revision: 1.1.10.3 $ */
///* The computational routine */
/*
to test the correlation routines


- uses the Wigner-d symmetries
- uses part of SpharmonicKit
- INTERLEAVED (i.e. real/imaginary) SAMPLES of signal and pattern files
- [result] -> optional -> filename of all the correlation values
(if you want all of them)
- bwIn -> bw of input spherical signals
- bwOut -> bw of so(3) transform you want to do
- degLim -> max degree of Wigner-D functions you'll be using


ASSUMES bwIn >= bwOut

example: test_soft_sym_correlate2 signalFile patternFile bwIn bwOut degLim [result]

*/

#include "mex.h"
#include <math.h>
#define M_PI       3.14159265358979323846

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "so3_correlate_sym.h"
#include "soft_sym.h"

#include "FST_semi_memo.h" 
#include "cospmls.h"
#include "s2_primitive.h"
#include "legendreTransforms.h"


void  ComputeCorrelationOnSphere(double* dSignalRealPart, double* dSignalImagPart, double* dPatternRealPart, double* dPatternImagPart, double bandwidthForward, double bandwidthInverse, double degrees, double* alpha, double* beta, double* gamma)
{
	// Function extracted from minimal changes that are required for the correlation to work from external data source such as MATLAB
	//FILE *fp;
	int i;
	double n, bwIn, bwOut, degLim;
	//double tstart, tstop;
	double *workspace1, *workspace2;
	double *sigR, *sigI;
	double *sigCoefR, *sigCoefI;
	double *patCoefR, *patCoefI;
	double *so3SigR, *so3SigI;
	double *so3CoefR, *so3CoefI;
	int tmp, maxloc, ii, jj, kk;
	double maxval;
	double *seminaive_naive_tablespace;
	double **seminaive_naive_table;


	bwIn = bandwidthForward;
	bwOut = bandwidthInverse;
	degLim = degrees;

	n = 2 * bwIn;

	sigR = (double *)calloc(n * n, sizeof(double));
	sigI = (double *)calloc(n * n, sizeof(double));
	so3SigR = (double *)malloc(sizeof(double) * (8 * bwOut*bwOut*bwOut));
	so3SigI = (double *)malloc(sizeof(double) * (8 * bwOut*bwOut*bwOut));
	workspace1 = (double *)malloc(sizeof(double) * (16 * bwOut*bwOut*bwOut));
	workspace2 = (double *)malloc(sizeof(double) * ((14 * bwIn*bwIn) + (48 * bwIn)));
	sigCoefR = (double *)malloc(sizeof(double) * bwIn * bwIn);
	sigCoefI = (double *)malloc(sizeof(double) * bwIn * bwIn);
	patCoefR = (double *)malloc(sizeof(double) * bwIn * bwIn);
	patCoefI = (double *)malloc(sizeof(double) * bwIn * bwIn);
	so3CoefR = (double *)malloc(sizeof(double) * ((4 * bwOut*bwOut*bwOut - bwOut) / 3));
	so3CoefI = (double *)malloc(sizeof(double) * ((4 * bwOut*bwOut*bwOut - bwOut) / 3));

	seminaive_naive_tablespace =
		(double *)malloc(sizeof(double) *
		(Reduced_Naive_TableSize(bwIn, bwIn) +
		Reduced_SpharmonicTableSize(bwIn, bwIn)));


	/****
	At this point, check to see if all the memory has been
	allocated. If it has not, there's no point in going further.
	****/

	if ((seminaive_naive_tablespace == NULL) ||
		(sigR == NULL) || (sigI == NULL) ||
		(so3CoefR == NULL) || (so3CoefI == NULL) ||
		(workspace1 == NULL) || (workspace2 == NULL) ||
		(sigCoefR == NULL) || (sigCoefI == NULL) ||
		(patCoefR == NULL) || (patCoefI == NULL) ||
		(so3CoefR == NULL) || (so3CoefI == NULL))
	{
		perror("Error in allocating memory");
		exit(1);
	}

	mexPrintf( "Generating seminaive_naive tables...\n");
	seminaive_naive_table = SemiNaive_Naive_Pml_Table(bwIn, bwIn,
		seminaive_naive_tablespace,
		workspace2);

	mexPrintf("now taking spherical transform of signal\n");
	FST_semi_memo(dSignalRealPart, dSignalImagPart,
		sigCoefR, sigCoefI,
		n, seminaive_naive_table,
		workspace2, 1, bwIn);
	if (1)
		for (i = 0; i<5; i++)
			mexPrintf("%d\t%f\t%f\n", i, sigCoefR[i], sigCoefI[i]);


	mexPrintf("now taking spherical transform of pattern\n");
	FST_semi_memo(dPatternRealPart, dPatternImagPart,
		patCoefR, patCoefI,
		n, seminaive_naive_table,
		workspace2, 1, bwIn);

	if (1)
		for (i = 0; i<5; i++)
			mexPrintf("%d\t%f\t%f\n", i, patCoefR[i], patCoefI[i]);

	mexPrintf("freeing seminaive_naive_table and seminaive_naive_tablespace\n");

	free(seminaive_naive_table);
	free(seminaive_naive_tablespace);


	mexPrintf("about to combine coefficients\n");

	/* combine coefficients */

	//clock_t start_t, end_t, total_t;


	//tstart = clock() / CLOCKS_PER_SEC; //  csecond();
	so3CombineCoef(bwIn, bwOut, degLim,
		sigCoefR, sigCoefI,
		patCoefR, patCoefI,
		so3CoefR, so3CoefI);
	//tstop = clock() / CLOCKS_PER_SEC;// csecond();
	//mexPrintf(stderr, "combine time \t = %.4e\n", tstop - tstart);

	mexPrintf("about to inverse so(3) transform\n");

	//tstart = clock() / CLOCKS_PER_SEC;// csecond();
	/* now inverse so(3) */
	Inverse_SO3_Naive_sym(bwOut,
		so3CoefR, so3CoefI,
		so3SigR, so3SigI,
		workspace1, workspace2,
		1);
	//tstop = clock() / CLOCKS_PER_SEC; //csecond();
	mexPrintf("finished inverse so(3) transform\n");
	//mexPrintf(stderr, "inverse so(3) time \t = %.4e\n", tstop - tstart);


	/* now find max value */
	maxval = 0.0;
	maxloc = 0;
	for (i = 0; i < 8 * bwOut*bwOut*bwOut; i++)
	{
		if (so3SigR[i] >= maxval)
		{
			maxval = so3SigR[i];
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

	
	free(so3CoefI);
	free(so3CoefR);
	free(patCoefI);
	free(patCoefR);
	free(sigCoefI);
	free(sigCoefR);
	free(workspace2);
	free(workspace1);

	free(so3SigI);
	free(so3SigR);

	free(sigI);
	free(sigR);
}


/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* check for proper number of arguments */
	if (nrhs != 7) {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "Seven inputs required.");
	}
	if (nlhs != 3) {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs", "Three output required.");
	}

	// Actual Processing
	double* dSignalRealPart;
	double*	dSignalImagPart;
	double* dPatternRealPart;
	double* dPatternImagPart;    // Input Arrays
	double bandwidthForward, bandwidthInverse, degrees;                                 // Input scalars
	double* alpha = 0;
	double* beta = 0;
	double* gamma =0;														// Output scalars		

	/* create a pointer to the real data in the input matrix  */
	dSignalRealPart = mxGetPr(prhs[0]);
	dSignalImagPart = mxGetPr(prhs[1]);
	dPatternRealPart = mxGetPr(prhs[2]);
	dPatternImagPart = mxGetPr(prhs[3]);

	/* get the value of the scalar input  */
	bandwidthForward = mxGetScalar(prhs[4]);
	bandwidthInverse = mxGetScalar(prhs[5]);
	degrees = mxGetScalar(prhs[6]);

	/* Assign a pointer to the output */
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);

	
	alpha = mxGetPr(plhs[0]);
	beta = mxGetPr(plhs[1]);
	gamma = mxGetPr(plhs[2]);

	/* call the computational routine */
	ComputeCorrelationOnSphere(dSignalRealPart, dSignalImagPart, dPatternRealPart, dPatternImagPart, bandwidthForward, bandwidthInverse, degrees, alpha, beta, gamma);

}
