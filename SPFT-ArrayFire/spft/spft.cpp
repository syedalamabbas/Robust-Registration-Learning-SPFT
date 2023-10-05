/*
* This is a program for transforming Volumetric Spatial Grid to Spherical Polar Fourier Transform Spectral Grid
*/

#include <arrayfire.h>
#include <af/util.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <windows.h>
#include <tlhelp32.h>
#include <stdio.h>

#include <filesystem>
namespace fs = std::filesystem;

using namespace af;

#include <ppl.h>
using namespace concurrency;

enum AXES
{
	X_AXIS = 1,
	Y_AXIS = 0,   // This is how arrayFire refers to Y-axis
	Z_AXIS = 2
};

static const af::dtype PRECISION_REAL = f64;                   // Change here from single to double
static const af::dtype PRECISION_COMPLEX = c64;

array cexp(const array &in)
{
	if (!in.iscomplex()) return exp(in);
	return exp(real(in))*complex(cos(imag(in)), sin(imag(in)));
}

void SplitMultiplyComplex(array& A, array&  B, array& realRealPart, array& realImagPart, array& imagRealPart, array& imagImagPart)
{
	// Consider multiplication of complex numbers A = (a+ib); B = (c+id)	
	//array A, B;
	//array realRealPart;             // ac
	//array realImagPart;				// ad
	//array imagRealPart;				// bc
	//array imagImagPart;				// bd
	realRealPart = real(A)*real(B);
	//af_print(realRealPart);
	realImagPart = real(A)*imag(B);
	imagRealPart = imag(A)*real(B);
	imagImagPart = imag(A)*imag(B);
}

static array indexedElements_Tiled_En;
static array indexedElements_Tiled_PreMulti;
static array indexedElements_Tiled_PostMulti;

void GlobalArraysComputeInitialize(int d_NoOfElements)
{
	/*-------------------------------------------Creating Index Cubes and Sequences----------------------------------------------------*/
	int N = d_NoOfElements - 1;
	array leftSideIndexes = array(seq(0, N, 1)).as(PRECISION_REAL);
	// af_print(leftSideIndexes);
	array rightSideIndexesOnes = -1 * array(seq(1, d_NoOfElements, 1)).as(PRECISION_REAL);
	//af_print(flip( rightSideIndexesOnes,0));
	array rightSideIndexesZeros = constant(0, d_NoOfElements, 1, PRECISION_REAL);
	//af_print(rightSideIndexesZeros);
	array rightSideIndexesN_2 = constant(N / 2, d_NoOfElements, 1, PRECISION_REAL);
	//af_print(rightSideIndexesN_2);

	array indexedElementsEn = join(0, leftSideIndexes, flip(rightSideIndexesOnes, 0));
	array indexedElementsPre = join(0, leftSideIndexes, rightSideIndexesN_2);     /* This is for Keeping pre and post multiplication factor upper half only*/
	array indexedElementsPost = join(0, leftSideIndexes, rightSideIndexesZeros);
	//af_print(indexedElementsEn);
	//af_print(indexedElementsPre_Post);

	 indexedElements_Tiled_En = tile(pow(indexedElementsEn, 2), 1, d_NoOfElements, d_NoOfElements);
	//af_print(indexedElements_Tiled_En.slice(0));
	 indexedElements_Tiled_PreMulti = tile(indexedElementsPre - N / 2, 1, d_NoOfElements, d_NoOfElements);
	//af_print(indexedElements_Tiled_PreMulti.slice(0));
	 indexedElements_Tiled_PostMulti = tile(indexedElementsPost, 1, d_NoOfElements, d_NoOfElements);
	//af_print(indexedElements_Tiled_PostMulti.slice(0));
}

array ComputeNonVectorColumnwise_FrFTUniform(array& Image_3D, double ColumnScale, int d_NoOfElements)
{
	// Image_3D is unpadded just a 3D volume (d_NoOfElements)^3 whose column wise FrFT scaling needs to be done
	// ColumnScale is the given scale that needs to be operated on the image uniformly for entire volume (3D) 
	// d_NoOfElements is the number of elements in the volume Image_3D which is a perfect cube, 


	/**************************************************Preparing the 3D Volume Stack***************************************************/
	array zeros = constant(0, d_NoOfElements, d_NoOfElements, d_NoOfElements, PRECISION_REAL);
	array paddedImage_3D = join(0, Image_3D, zeros);									  // Column Wise padding
	//af_print(paddedImage_3D.slice(0));
	int N = d_NoOfElements - 1;

	/*-------------------------------------------Creating Index Cubes and Sequences----------------------------------------------------*/
	
	//array leftSideIndexes = array(seq(0, N, 1)).as(PRECISION_REAL);
	//// af_print(leftSideIndexes);
	//array rightSideIndexesOnes = -1 * array(seq(1, d_NoOfElements, 1)).as(PRECISION_REAL);
	////af_print(flip( rightSideIndexesOnes,0));
	//array rightSideIndexesZeros = constant(0, d_NoOfElements, 1, PRECISION_REAL);
	////af_print(rightSideIndexesZeros);
	//array rightSideIndexesN_2 = constant(N / 2, d_NoOfElements, 1, PRECISION_REAL);
	////af_print(rightSideIndexesN_2);

	//array indexedElementsEn = join(0, leftSideIndexes, flip(rightSideIndexesOnes, 0));
	//array indexedElementsPre = join(0, leftSideIndexes, rightSideIndexesN_2);     /* This is for Keeping pre and post multiplication factor upper half only*/
	//array indexedElementsPost = join(0, leftSideIndexes, rightSideIndexesZeros);
	////af_print(indexedElementsEn);
	////af_print(indexedElementsPre_Post);

	//array indexedElements_Tiled_En = tile(pow(indexedElementsEn, 2), 1, d_NoOfElements, d_NoOfElements);
	////af_print(indexedElements_Tiled_En.slice(0));
	//array indexedElements_Tiled_PreMulti = tile(indexedElementsPre - N / 2, 1, d_NoOfElements, d_NoOfElements);
	////af_print(indexedElements_Tiled_PreMulti.slice(0));
	//array indexedElements_Tiled_PostMulti = tile(indexedElementsPost, 1, d_NoOfElements, d_NoOfElements);
	////af_print(indexedElements_Tiled_PostMulti.slice(0));


	/*-------------------Precomputing the Essential Sequence Cubes :: All complex values here --------------------*/
	//array ImaginaryUnitTiled = tile(i, 2 * d_NoOfElements, d_NoOfElements, d_NoOfElements);
	//af_print(ImaginaryUnitTiled.slice(0));
	array Seq_En = cexp(complex(0, -Pi * indexedElements_Tiled_En * ColumnScale / d_NoOfElements));   /* E(n) as defined in the paper*/
	//af_print((Pi * indexedElements_Tiled_En * ColumnScale / d_NoOfElements).slice(0));
	//af_print(Seq_En.slice(0));
	//array ones = constant(1, d_NoOfElements, 1, PRECISION_REAL);
	//array subtractValues = tile(join(0, rightSideIndexesZeros, ones), 1, d_NoOfElements, d_NoOfElements);			/* This is for  Keeping pre and post multiplication factor upper half only*/
	array PreMultiplicationFactor = cexp(complex(0, af::Pi * indexedElements_Tiled_PreMulti * ColumnScale * N / d_NoOfElements));// -subtractValues;
	//af_print(PreMultiplicationFactor.slice(0));
	array PostMultiplicationFactor = cexp(complex(0, af::Pi *  indexedElements_Tiled_PostMulti * ColumnScale * N / d_NoOfElements));// - subtractValues;
	//af_print(PostMultiplicationFactor.slice(0));

	/*--------------------Preprocessing Cubes-----------------------*/
	array Image_PreMulti = paddedImage_3D * PreMultiplicationFactor;
	//af_print(Image_PreMulti.slice(0));
	array Image_PreMulti_SeqEn = Image_PreMulti * Seq_En;
	//af_print(Image_PreMulti_SeqEn.slice(0));

	/*-------------------Computing Convolution--------------------*/
	array firstFFT_Image = fft(Image_PreMulti_SeqEn);
	//af_print(firstFFT_Image.slice(0));
	array secondFFT_Seq_En = fft(conjg(Seq_En));
	//af_print(secondFFT_Seq_En.slice(0));
	array interim_FrFT_Image = ifft(firstFFT_Image * secondFFT_Seq_En);
	//af_print(interim_FrFT_Image.slice(0));

	/*-------------------Postprocessing-----------------------------*/
	array  Final_FrFT_Image = interim_FrFT_Image * Seq_En * PostMultiplicationFactor;
	//af_print(Final_FrFT_Image.slice(0));

	/*--------------------Grab only the top half elements drop overlapping------------------*/
	Final_FrFT_Image = Final_FrFT_Image.rows(0, N);
	//af_print(Final_FrFT_Image.slice(0));

	return Final_FrFT_Image;
}

void ComputeNonVector_NonUniformFrFT_SingleColumn(array&FrFT1D_Image_3D, double ColumnScale, int d_NoOfElements, array& Slice2D, array& Slice2D_Conj)
{
	int N = d_NoOfElements - 1;
	array lineSpacing = array(seq(-N / 2, N / 2)).as(PRECISION_REAL);
	array lineSpacing_tiled_Y = tile(lineSpacing, 1, d_NoOfElements);
	//af_print(lineSpacing_tiled_Y);
	array lineSpacing_Square = lineSpacing_tiled_Y * lineSpacing_tiled_Y.T();
	//af_print(lineSpacing_Square);
	array lineSpacing_Square_TiledLevel = tile(lineSpacing_Square, 1, 1, d_NoOfElements);
	//af_print(lineSpacing_Square_TiledLevel);
	array Column_ScaleFactor = cexp(complex(0, -2 * af::Pi * lineSpacing_Square_TiledLevel * ColumnScale / d_NoOfElements));
	//af_print(Column_ScaleFactor.slice(0))

	// Consider multiplication of complex numbers A = (a+ib); B = (c+id)	
	array realRealPart;             // ac
	array realImagPart;				// ad 
	array imagRealPart;				// bc
	array imagImagPart;				// bd

	SplitMultiplyComplex(FrFT1D_Image_3D, Column_ScaleFactor, realRealPart, realImagPart, imagRealPart, imagImagPart);
	Slice2D = sum(complex(realRealPart - imagImagPart, realImagPart + imagRealPart)); // sum(flip(FrFT_Image_X_Cube, 0) *BetaFactor);
	Slice2D_Conj = sum(complex(realRealPart + imagImagPart, imagRealPart - realImagPart));//  sum(flip(FrFT_Image_X_Cube, 0) * conjg(BetaFactor));

	/*af_print(Slice2D);
	af_print(Slice2D_Conj);*/

	Slice2D = moddims(Slice2D, d_NoOfElements, d_NoOfElements);                // Making it 2D without changing the order of data
	Slice2D_Conj = moddims(Slice2D_Conj, d_NoOfElements, d_NoOfElements);

	/*af_print(Slice2D);
	af_print(Slice2D_Conj);*/
}

void ComputeFinal_2ComplementaryLines(array& FrFTFrFT_Image2D, double ColumnScale, int d_NoOfElements, array& Line1D, array& Line1D_Conj)
{
	int N = d_NoOfElements - 1;
	array lineSpacing = array(seq(-N / 2, N / 2)).as(PRECISION_REAL);
	array lineSpacing_tiled_Y = tile(lineSpacing, 1, d_NoOfElements);
	//af_print(lineSpacing_tiled_Y);
	array lineSpacing_Square = lineSpacing_tiled_Y * lineSpacing_tiled_Y.T();
	array Column_ScaleFactor = cexp(complex(0, -2 * af::Pi * lineSpacing_Square * ColumnScale / d_NoOfElements));
	//af_print(Column_ScaleFactor)

	// Consider multiplication of complex numbers A = (a+ib); B = (c+id)	
	array realRealPart;             // ac
	array realImagPart;				// ad 
	array imagRealPart;				// bc
	array imagImagPart;				// bd

	SplitMultiplyComplex(FrFTFrFT_Image2D, Column_ScaleFactor, realRealPart, realImagPart, imagRealPart, imagImagPart);
	Line1D = sum(complex(realRealPart - imagImagPart, realImagPart + imagRealPart)); // sum(flip(FrFT_Image_X_Cube, 0) *BetaFactor);
	Line1D_Conj = sum(complex(realRealPart + imagImagPart, imagRealPart - realImagPart));//  sum(flip(FrFT_Image_X_Cube, 0) * conjg(BetaFactor));
	//af_print(Line1D);
	//af_print(Line1D_Conj);
}

array ComputeForSinglePoint3D(array& Image_3D, array& gridSpacing, int d_NoOfElements, double* desiredPoint)
{
	int N = d_NoOfElements - 1;

	// First since Columns are first dimension in ArrayFire , Y-axis scaling first
	array y_Multiplicationfactor = cexp(complex(0, -2 * Pi*desiredPoint[1] * gridSpacing / d_NoOfElements));

	array Map_y = tile(y_Multiplicationfactor, 1, d_NoOfElements, d_NoOfElements);
	//af_print(Image_3D);
	//af_print(Map_y);
	array Image_2D = moddims(sum(Image_3D*Map_y), d_NoOfElements, d_NoOfElements);
	//af_print(Image_2D);

	// X-axis 

	array x_Multiplicationfactor = cexp(-complex(0, 2 * Pi*desiredPoint[0] * gridSpacing / d_NoOfElements));

	array Map_x = tile(x_Multiplicationfactor, 1, d_NoOfElements);
	//af_print(Map_x);
	array Image_1D = sum(Image_2D*Map_x);
	//af_print(Image_1D);

	// Finally Z-axis 

	array z_Multiplicationfactor = cexp(-complex(0, 2 * Pi*desiredPoint[2] * gridSpacing / d_NoOfElements));
	Image_1D = moddims(Image_1D, d_NoOfElements, 1);
	//af_print(z_Multiplicationfactor);
	//af_print(Image_1D);
	array SinglePoint = sum(Image_1D* z_Multiplicationfactor);
	return SinglePoint;

}

array Get2DFullPolarDFT(array& Image, int d_NoOfElements, int d_NoOfAngles, int d_NoOfLevels)    // This will require the same processing as defined for 2D see the solution in 2D
{
	
	int N = d_NoOfElements - 1;
	cdouble i_cdouble = { 0, 1 };
	array i = constant(i_cdouble, 1, 1, PRECISION_COMPLEX);/* imaginary unit */

	/*-----------------------------------Preparing Padded & Tiled Images in X and Y--------------------------------------------*/
	//array Zeros = constant(0, d_NoOfElements, d_NoOfElements, PRECISION_REAL);
	array Zeros = constant(0, d_NoOfElements, d_NoOfElements, PRECISION_COMPLEX);
	array Ones = constant(1, d_NoOfElements, d_NoOfElements, PRECISION_COMPLEX);
	//af_print(Image);
	//af_print(Zeros);
	array Zero_Padded_X_Image_Transposed = join(1, Image, Zeros).T();        /* Transposing so that the rows acts as columns to be operated upon for FFTs*/
	array Zero_Padded_Y_Image = join(0, flip(Image, 0), Zeros);              /* Needed since image gathered gives +ve to -ve indexes and we require opposite for computations*/
	array Image_Tiled_X = tile(Zero_Padded_X_Image_Transposed, 1, 1, d_NoOfLevels);
	array Image_Tiled_Y = tile(Zero_Padded_Y_Image, 1, 1, d_NoOfLevels);

	/*-------------------------------------------Creating Index Cubes and Sequences----------------------------------------------------*/
	array leftSideIndexes = array(seq(0, N, 1)).as(PRECISION_REAL);
	array rightSideIndexesOnes = -1 * array(seq(1, d_NoOfElements, 1)).as(PRECISION_REAL);
	array rightSideIndexesZeros = constant(0, d_NoOfElements, 1, PRECISION_REAL);
	array rightSideIndexesN_2 = constant(N / 2, d_NoOfElements, 1, PRECISION_REAL);

	array indexedElementsEn = join(0, leftSideIndexes, flip(rightSideIndexesOnes, 0));
	array indexedElementsPre = join(0, leftSideIndexes, rightSideIndexesN_2);     /* This is for Keeping pre and post multiplication factor upper half only*/
	array indexedElementsPost = join(0, leftSideIndexes, rightSideIndexesZeros);

	array indexedElements_Tiled_En = tile(pow(indexedElementsEn, 2), 1, d_NoOfElements, d_NoOfLevels);
	array indexedElements_Tiled_PreMulti = tile(indexedElementsPre - N / 2, 1, d_NoOfElements, d_NoOfLevels);
	array indexedElements_Tiled_PostMulti = tile(indexedElementsPost, 1, d_NoOfElements, d_NoOfLevels);

	array lineSpacing = array(seq(-N / 2, N / 2)).as(PRECISION_REAL);
	array lineSpacing_tiled_Y = tile(lineSpacing, 1, d_NoOfElements);
	array lineSpacing_Square = lineSpacing_tiled_Y * lineSpacing_tiled_Y.T();
	array lineSpacing_Square_TiledLevel = tile(lineSpacing_Square, 1, 1, d_NoOfLevels);
	
	/*--------------------------Creating Alpha Cubes------------------------------------*/
	array alpha_Levels = cos(af::Pi / d_NoOfAngles * array(seq(1, d_NoOfLevels)).as(PRECISION_REAL));
	array alpha_Mods = moddims(alpha_Levels, 1, 1, d_NoOfLevels);
	array alpha_Tiled_depth = tile(alpha_Mods, 2 * d_NoOfElements, d_NoOfElements, 1);
	array alpha_Tiled = moddims(alpha_Tiled_depth, 2 * d_NoOfElements, d_NoOfElements, d_NoOfLevels);


	/*--------------------------Creating Beta Cubes---------------------------------------*/
	array beta_Levels = sin(Pi / d_NoOfAngles * array(seq(1, d_NoOfLevels)).as(PRECISION_REAL));
	array beta_Mods = moddims(beta_Levels, 1, 1, d_NoOfLevels);
	array beta_Tiled_depth = tile(beta_Mods, d_NoOfElements, d_NoOfElements, 1);
	array beta_Tiled = moddims(beta_Tiled_depth, d_NoOfElements, d_NoOfElements, d_NoOfLevels);
	

	/*-------------------Precomputing the Essential Sequence Cubes :: All complex values here --------------------*/
	array Seq_En = cexp(-tile(i, 2 * d_NoOfElements, d_NoOfElements, d_NoOfLevels) * af::Pi * indexedElements_Tiled_En * alpha_Tiled / d_NoOfElements);   /* E(n) as defined in the paper*/
	
	array subtractValues = tile(join(0, Zeros, Ones), 1, 1, d_NoOfLevels);			/* This is for  Keeping pre and post multiplication factor upper half only*/
	array PreMultiplicationFactor = cexp(tile(i, 2 * d_NoOfElements, d_NoOfElements, d_NoOfLevels) * af::Pi * indexedElements_Tiled_PreMulti * alpha_Tiled * N / d_NoOfElements) - subtractValues;
	array PostMultiplicationFactor = cexp(tile(i, 2 * d_NoOfElements, d_NoOfElements, d_NoOfLevels)* af::Pi *  indexedElements_Tiled_PostMulti * alpha_Tiled * N / d_NoOfElements) - subtractValues;
	array BetaFactor = cexp(-2 * tile(i, d_NoOfElements, d_NoOfElements, d_NoOfLevels) * af::Pi * lineSpacing_Square_TiledLevel * beta_Tiled / d_NoOfElements);

	array ZeroNinty_Factor = cexp(-2 * tile(i, d_NoOfElements, d_NoOfElements) * af::Pi * lineSpacing_Square * 1 / d_NoOfElements);

	/*--------------------Preprocessing Cubes-----------------------*/
	array Image_Tiled_X_PreMulti = Image_Tiled_X * PreMultiplicationFactor;
	array Image_Tiled_X_PreMulti_SeqEn = Image_Tiled_X_PreMulti * Seq_En;
	array Image_Tiled_Y_PreMulti = Image_Tiled_Y * PreMultiplicationFactor;
	array Image_Tiled_Y_PreMulti_SeqEn = Image_Tiled_Y_PreMulti * Seq_En;


	/*-------------------Computing Convolution--------------------*/
	int normalizationFactor = pow(2 * d_NoOfElements, 0);			// This is a scalar single value especially needed only for ArrayFire ifft
	array firstFFT_X = fft(moddims(Image_Tiled_X_PreMulti_SeqEn, 2 * d_NoOfElements, d_NoOfLevels*d_NoOfElements));
	array secondFFT_X = fft(moddims(conjg(Seq_En), 2 * d_NoOfElements, d_NoOfLevels*d_NoOfElements));
	array interim_FrFT_X = ifft(firstFFT_X * secondFFT_X) / normalizationFactor;
	array firstFFT_Y = fft(moddims(Image_Tiled_Y_PreMulti_SeqEn, 2 * d_NoOfElements, d_NoOfLevels*d_NoOfElements));
	array secondFFT_Y = secondFFT_X;
	array interim_FrFT_Y = ifft(firstFFT_Y * secondFFT_Y) / normalizationFactor;


	interim_FrFT_X = moddims(interim_FrFT_X, 2 * d_NoOfElements, d_NoOfElements, d_NoOfLevels);    /* Rearranging matrix to level based volume*/
	interim_FrFT_Y = moddims(interim_FrFT_Y, 2 * d_NoOfElements, d_NoOfElements, d_NoOfLevels);

	/*-------------------Postprocessing-----------------------------*/
	array  FrFT_Image_X = interim_FrFT_X * Seq_En * PostMultiplicationFactor;
	array  FrFT_Image_Y = interim_FrFT_Y * Seq_En * PostMultiplicationFactor;


	/*--------------------Grab only the top half elements drop overlapping------------------*/
	array FrFT_Image_X_Cube = FrFT_Image_X.rows(0, N);
	array FrFT_Image_Y_Cube = FrFT_Image_Y.rows(0, N);

	/*--------------------Finally all computations for  the  Polar Grid-----------*/
	//   Computing for all the grid expect two special indexes
	array levelSeq = array(seq(0, d_NoOfLevels - 1)).as(PRECISION_REAL);
	array finalIndexSeq1_X = 1 + levelSeq;
	array finalIndexSeq2_X = d_NoOfAngles - finalIndexSeq1_X;
	array finalIndexSeq3_Y = d_NoOfAngles / 2 - finalIndexSeq1_X;
	array finalIndexSeq4_Y = d_NoOfAngles / 2 + finalIndexSeq1_X;

	/* This transpose is Required since we will no operate in the column wise transposed axis */
	FrFT_Image_X_Cube = FrFT_Image_X_Cube.T();
	FrFT_Image_Y_Cube = FrFT_Image_Y_Cube.T();

	// Consider multiplication of complex numbers A = (a+ib); B = (c+id), This is used to reduce the computations by half 	
	array realRealPart;             // ac
	array realImagPart;				// ad
	array imagRealPart;				// bc
	array imagImagPart;				// bd

	SplitMultiplyComplex(flip(FrFT_Image_X_Cube, 0), BetaFactor, realRealPart, realImagPart, imagRealPart, imagImagPart);
	array tempSeq_X = sum(complex(realRealPart - imagImagPart, realImagPart + imagRealPart)); // sum(flip(FrFT_Image_X_Cube, 0) *BetaFactor);
	array tempSeqConj_X = sum(complex(realRealPart + imagImagPart, imagRealPart - realImagPart));//  sum(flip(FrFT_Image_X_Cube, 0) * conjg(BetaFactor));
	SplitMultiplyComplex(FrFT_Image_Y_Cube, BetaFactor, realRealPart, realImagPart, imagRealPart, imagImagPart);
	array tempSeq_Y = sum(complex(realRealPart - imagImagPart, realImagPart + imagRealPart)); //sum(FrFT_Image_Y_Cube * BetaFactor);
	array tempSeqConj_Y = sum(complex(realRealPart + imagImagPart, imagRealPart - realImagPart)); //sum(FrFT_Image_Y_Cube *conjg(BetaFactor));

	array finalSeq_X = moddims(tempSeq_X, d_NoOfElements, d_NoOfLevels).T();
	array finalSeqConj_X = moddims(tempSeqConj_X, d_NoOfElements, d_NoOfLevels).T();
	finalSeqConj_X = flip(finalSeqConj_X, 1);
	array finalSeq_Y = moddims(tempSeq_Y, d_NoOfElements, d_NoOfLevels).T();
	array finalSeqConj_Y = moddims(tempSeqConj_Y, d_NoOfElements, d_NoOfLevels).T();

	// Removing just 2 redundant computations for 45 degree case
	if (0 == remainder(d_NoOfAngles, 4))
	{
		finalIndexSeq3_Y = finalIndexSeq3_Y.rows(0, d_NoOfLevels - 2);          // Removing just the last rows from 4 structures
		finalSeq_Y = finalSeq_Y.rows(0, d_NoOfLevels - 2);
		finalIndexSeq4_Y = finalIndexSeq4_Y.rows(0, d_NoOfLevels - 2);
		finalSeqConj_Y = finalSeqConj_Y.rows(0, d_NoOfLevels - 2);
	}

	//   Computing seperately for two special indexes
	array SpecialTwoIndexes;
	array ZeroLineFrFT_Image_X_Cube;
	array NintyLineFrFT_Image_Y_Cube;
	if (PRECISION_REAL == f32)
	{
	   float zeroIndex = 0;
	   float nintyIndex = d_NoOfAngles / 2;
	    float values[] = { zeroIndex, nintyIndex };
	   SpecialTwoIndexes= array(2, 1, values);
	   ZeroLineFrFT_Image_X_Cube = FrFT_Image_Y_Cube.slice(zeroIndex).col(N / 2);
	   NintyLineFrFT_Image_Y_Cube = FrFT_Image_X_Cube.slice(zeroIndex).col(N / 2);
	}
	else
	{
		double zeroIndex = 0;
		double nintyIndex = d_NoOfAngles / 2;
		double values[] = { zeroIndex, nintyIndex };
		SpecialTwoIndexes = array(2, 1, values);
		ZeroLineFrFT_Image_X_Cube = FrFT_Image_Y_Cube.slice(zeroIndex).col(N / 2);
		NintyLineFrFT_Image_Y_Cube = FrFT_Image_X_Cube.slice(zeroIndex).col(N / 2);
	}

	

	array DFTZeroLine = sum(tile(ZeroLineFrFT_Image_X_Cube, 1, d_NoOfElements) *ZeroNinty_Factor);
	array DFTNinetyLine = sum(tile(flip(NintyLineFrFT_Image_Y_Cube, 0), 1, d_NoOfElements) *ZeroNinty_Factor);
	array SpecialTwoLines = join(0, DFTZeroLine, DFTNinetyLine);


	array UnsortedIndexes = join(0, join(0, join(0, join(0, finalIndexSeq1_X, finalIndexSeq2_X), finalIndexSeq3_Y), finalIndexSeq4_Y), SpecialTwoIndexes);
	array tiledUnsortedIndexes = tile(UnsortedIndexes, 1, d_NoOfElements);
	array UnsortedPolarGrid = join(0, join(0, join(0, join(0, finalSeq_X, finalSeqConj_X), finalSeq_Y), finalSeqConj_Y), SpecialTwoLines);


	array FinalPolarGridReal;// = constant(0, d_NoOfElements, d_NoOfAngles, PRECISION_COMPLEX);
	array Output_Keys_Sorted;
	sort(Output_Keys_Sorted, FinalPolarGridReal, tiledUnsortedIndexes, real(UnsortedPolarGrid));

	array FinalPolarGridImag;// = constant(0, d_NoOfElements, d_NoOfAngles, PRECISION_COMPLEX);
	array Output_Keys_Sorted2;
	sort(Output_Keys_Sorted2, FinalPolarGridImag, tiledUnsortedIndexes, imag(UnsortedPolarGrid));


	return complex(FinalPolarGridReal, FinalPolarGridImag);
}

static void ComputeFastSphericalPolarFourierTransform(array& Image3D, array& FinalSphericalGrid, int d_NoOfElements, int d_NoOfAnglesTheta, int d_NoOfLevelsTheta, int d_NoOfAnglesPhi, int d_NoOfLevelsPhi)
{
	array ReorderedImage_OperateColumns;                 // Reordered 3D Image at the start of the FrFT operation  
	array FrFT1D_Uniform_Image3D;                        // 3D Image obtained after passing 1st stage of FrFT, all columns have now been operated on
	array Silce2D, Slice2D_Conj;                         // 2D Slices obtained at the 2nd stage of FrFT
	array Line1, Line1_Conj, Line2, Line2_Conj;			 // 1D Lines at the end of the 3rd stage of FrFT, which are final	
	double alpha_factor, beta_factor, gamma_factor;      // Scaling factors in X-axis , Y-axis and Z-axis respectively, which change depending on the computation of the block 
	
	int N = d_NoOfElements - 1;
	GlobalArraysComputeInitialize(d_NoOfElements);
	for (int q = 1; q <= d_NoOfLevelsPhi; ++q)                // Splitting as defined in the paper (See the schematic of spherical divisions)
	{
		double anglePhi = q * Pi / d_NoOfAnglesPhi;

		//timer::start();
		for (int p = 1; p <= d_NoOfLevelsTheta; ++p)       // X oriented pair of Polar slices, ( q, d_NoOfAnglesPhi - q)
		{
			double angleTheta = p  * Pi / d_NoOfAnglesTheta;

			/**************************************************************************************************************************************************************************/
			// XX block  -- Concentric rectangles in YZ tiled along X -axis
			alpha_factor = cos(angleTheta) * cos(anglePhi);											         // Scaling needed as defined in the paper
			beta_factor = cos(angleTheta) * sin(anglePhi);
			gamma_factor = sin(angleTheta);

			ReorderedImage_OperateColumns = reorder(Image3D, X_AXIS, Y_AXIS, Z_AXIS);                                //  1st Swap, X-axis as columns now
			//af_print(Image3D.T().slices(0, 2));
			//af_print(ReorderedImage_OperateColumns.slices(0,2));
			//ReorderedImage_OperateColumns = reorder(ReorderedImage_OperateColumns, X_AXIS, Y_AXIS, Z_AXIS );            // Restoring the image as it was
			//af_print(ReorderedImage_OperateColumns.slices(0,2));

			FrFT1D_Uniform_Image3D = ComputeNonVectorColumnwise_FrFTUniform(ReorderedImage_OperateColumns, alpha_factor, d_NoOfElements); // This is the most expensive operation of order (N+1)^3 log(N+1)
			//af_print(FrFT1D_Uniform_Image3D)
			if (q == 1 && p == 1)         // Computing Polar Slice at anglephi = 90 , only need to be computed once
			{
				array Central2D_YZSlice2D = FrFT1D_Uniform_Image3D(N / 2, span, span);
				Central2D_YZSlice2D = moddims(Central2D_YZSlice2D, d_NoOfElements, d_NoOfElements);
				//af_print(Central2D_YZSlice2D);

				array properOrientedSlice = flip(Central2D_YZSlice2D.T(), 0);          // Verified match 
				//af_print(properOrientedSlice);
				array Polar2D = Get2DFullPolarDFT(properOrientedSlice, d_NoOfElements, d_NoOfAnglesTheta, d_NoOfLevelsTheta);
				//af_print(FinalSphericalGrid(span, d_NoOfAnglesPhi / 2, span));
				Polar2D = moddims(Polar2D, d_NoOfAnglesTheta, 1, d_NoOfElements);
				//af_print(Polar2D);
				FinalSphericalGrid(span, d_NoOfAnglesPhi / 2, span) = Polar2D;
				//af_print(BruteForceSphericalGrid(span, d_NoOfAnglesPhi / 2, span) - FinalSphericalGrid(span, d_NoOfAnglesPhi / 2, span));
			}

			//af_print(FrFT1D_Uniform_Image3D.slice(0))
			//af_print((reorder(FrFT1D_Uniform_Image3D, Z_AXIS, Y_AXIS, X_AXIS)).slice(0))
			//af_print(FrFT1D_Uniform_Image3D);
			//af_print(reorder(FrFT1D_Uniform_Image3D, Y_AXIS, X_AXIS, Z_AXIS));
			ComputeNonVector_NonUniformFrFT_SingleColumn(reorder(FrFT1D_Uniform_Image3D, X_AXIS, Y_AXIS, Z_AXIS), beta_factor, d_NoOfElements, Silce2D, Slice2D_Conj); //  2nd Swap, Y-axis as columns
			//af_print(Silce2D);
			//af_print(Slice2D_Conj);

			// Gather 4 lines
			ComputeFinal_2ComplementaryLines(Silce2D.T(), gamma_factor, d_NoOfElements, Line1, Line1_Conj);          // 3rd Swap , Z axis as columns
			ComputeFinal_2ComplementaryLines(Slice2D_Conj.T(), gamma_factor, d_NoOfElements, Line2, Line2_Conj);      // 3rd Swap , Z axis as columns

			//af_print(Line1);
			//af_print(Line1_Conj); 
			//af_print(Line2);
			//af_print(Line2_Conj);

			FinalSphericalGrid(p, q, span) = Line1;
			FinalSphericalGrid(p, d_NoOfAnglesPhi - q, span) = conjg(Line2_Conj);
			FinalSphericalGrid(d_NoOfAnglesTheta - p, q, span) = conjg(Line1_Conj);											// 1 More Swap to match !
			FinalSphericalGrid(d_NoOfAnglesTheta - p, d_NoOfAnglesPhi - q, span) = Line2;

			/*printf("\n My Special Solution 1st XX block\n");
			af_print(moddims(FinalSphericalGrid(p, q, span), 1, d_NoOfElements));
			printf("\n");
			af_print(moddims(FinalSphericalGrid(p, d_NoOfAnglesPhi - q, span), 1, d_NoOfElements));
			printf("\n");
			af_print(moddims(FinalSphericalGrid(d_NoOfAnglesTheta - p, q, span), 1, d_NoOfElements));
			printf("\n");
			af_print(moddims(FinalSphericalGrid(d_NoOfAnglesTheta - p, d_NoOfAnglesPhi - q, span), 1, d_NoOfElements));
			printf("\n");*/
			/**************************************************************************************************************************************************************************/

			if (p * 180 / d_NoOfAnglesTheta == 45)       // At 45 degrees its redundant to compute in the next block so skip it
				continue;


			/**************************************************************************************************************************************************************************/
			// XZ block  -- Concentric rectangles in XY tiled along Z -axis
			alpha_factor = sin(angleTheta) * cos(anglePhi);											         // Scaling needed as defined in the paper
			beta_factor = sin(angleTheta) * sin(anglePhi);
			gamma_factor = cos(angleTheta);

			ReorderedImage_OperateColumns = reorder(Image3D, Z_AXIS, Y_AXIS, X_AXIS);                                   //  1st Swap, Z-axis elements in columns now
			//af_print(ReorderedImage_OperateColumns.slices(0,2));
			//ReorderedImage_OperateColumns = reorder(ReorderedImage_OperateColumns, X_AXIS, Z_AXIS, Y_AXIS );           // Restoring the image as it was
			//af_print(ReorderedImage_OperateColumns.slices(0, 2));
			FrFT1D_Uniform_Image3D = ComputeNonVectorColumnwise_FrFTUniform(ReorderedImage_OperateColumns, gamma_factor, d_NoOfElements); // This is the most expensive operation of order (N+1)^3 log(N+1)
			//af_print(FrFT1D_Uniform_Image3D);

			if (q == 1 && p == 1)         // Computing One special Polar Slice at angletheta = 0 , only need to be computed once
			{
				array Central2D_XYSlice2D = FrFT1D_Uniform_Image3D(N / 2, span, span);
				Central2D_XYSlice2D = moddims(Central2D_XYSlice2D, d_NoOfElements, d_NoOfElements);
				//af_print(Central2D_XYSlice2D);

				array properOrientedSlice = flip(Central2D_XYSlice2D, 0);          // Verified match 
				array Polar2D = Get2DFullPolarDFT(properOrientedSlice, d_NoOfElements, d_NoOfAnglesPhi, d_NoOfLevelsPhi);
				Polar2D = moddims(Polar2D, 1, d_NoOfAnglesPhi, d_NoOfElements);
				//af_print(Polar2D);
				FinalSphericalGrid(0, span, span) = Polar2D;
				//af_print(BruteForceSphericalGrid(0, span, span) - FinalSphericalGrid(0, span, span));
			}


			ComputeNonVector_NonUniformFrFT_SingleColumn(FrFT1D_Uniform_Image3D.T(), beta_factor, d_NoOfElements, Silce2D, Slice2D_Conj); //  2nd Swap, Y-axis as columns

			// Gather 4 lines
			ComputeFinal_2ComplementaryLines(Silce2D.T(), alpha_factor, d_NoOfElements, Line1, Line1_Conj);          // 3rd Swap , X axis as columns
			ComputeFinal_2ComplementaryLines(Slice2D_Conj.T(), alpha_factor, d_NoOfElements, Line2, Line2_Conj);     // Swap Lines ! // 3rd Swap , Z axis as columns

			FinalSphericalGrid(d_NoOfAnglesTheta / 2 - p, q, span) = Line1;
			FinalSphericalGrid(d_NoOfAnglesTheta / 2 - p, d_NoOfAnglesPhi - q, span) = (Line1_Conj);

			FinalSphericalGrid(d_NoOfAnglesTheta / 2 + p, q, span) = (Line2_Conj);
			FinalSphericalGrid(d_NoOfAnglesTheta / 2 + p, d_NoOfAnglesPhi - q, span) = (Line2);

			/*printf("\n My Special Solution 1st XZ block\n");
			af_print(moddims(FinalSphericalGrid(d_NoOfAnglesTheta / 2 - p, q, span), 1, d_NoOfElements));
			printf("\n");
			af_print(moddims(FinalSphericalGrid(d_NoOfAnglesTheta / 2 - p, d_NoOfAnglesPhi - q, span), 1, d_NoOfElements));
			printf("\n");
			af_print(moddims(FinalSphericalGrid(d_NoOfAnglesTheta / 2 + p, q, span), 1, d_NoOfElements));
			printf("\n");
			af_print(moddims(FinalSphericalGrid(d_NoOfAnglesTheta / 2 + p, d_NoOfAnglesPhi - q, span), 1, d_NoOfElements));
			printf("\n");*/

			/**************************************************************************************************************************************************************************/

			if (q * 180 / d_NoOfAnglesPhi == 45)       // At 45 degrees its redundant to compute same in the next block  so skip it
				continue;

			/**************************************************************************************************************************************************************************/
			// YZ block                                                   // Y oriented pair of Polar slices, (d_NoOfAnglesPhi / 2 - q, d_NoOfAnglesPhi / 2 + q)
			alpha_factor = sin(angleTheta) * sin(anglePhi);				 // Scaling needed as defined in the paper
			beta_factor = sin(angleTheta) * cos(anglePhi);
			gamma_factor = cos(angleTheta);

			// Reusing the previously computed FrFT 1D Uniform scaling since the gamma_factor for both XZ and YZ block is the same

			//ReorderedImage_OperateColumns = reorder(Image3D, Z_AXIS, Y_AXIS, X_AXIS);                                   //  1st Swap, Z-axis elements in columns now
			//FrFT1D_Uniform_Image3D = ComputeNonVectorColumnwise_FrFTUniform(ReorderedImage_OperateColumns, gamma_factor, d_NoOfElements); // This is the most expensive operation of order (N+1)/4*(N+1)^3 log(N+1)
			ComputeNonVector_NonUniformFrFT_SingleColumn(FrFT1D_Uniform_Image3D.T(), beta_factor, d_NoOfElements, Silce2D, Slice2D_Conj); //  2nd Swap, Y-axis as columns

			// Gather 4 lines
			ComputeFinal_2ComplementaryLines(Silce2D.T(), alpha_factor, d_NoOfElements, Line1, Line1_Conj);          // 3rd Swap , X axis as columns
			ComputeFinal_2ComplementaryLines(Slice2D_Conj.T(), alpha_factor, d_NoOfElements, Line2, Line2_Conj);     // Swap Lines ! // 3rd Swap , Z axis as columns

			FinalSphericalGrid(d_NoOfAnglesTheta / 2 - p, d_NoOfAnglesPhi / 2 - q, span) = Line1;
			FinalSphericalGrid(d_NoOfAnglesTheta / 2 - p, d_NoOfAnglesPhi / 2 + q, span) = Line1_Conj;

			FinalSphericalGrid(d_NoOfAnglesTheta / 2 + p, d_NoOfAnglesPhi / 2 - q, span) = Line2_Conj;
			FinalSphericalGrid(d_NoOfAnglesTheta / 2 + p, d_NoOfAnglesPhi / 2 + q, span) = Line2;

			/*printf("\n My Special Solution 1st YZ block\n");
			af_print(moddims(FinalSphericalGrid(d_NoOfAnglesTheta / 2 - p, d_NoOfAnglesPhi / 2 - q, span), 1, d_NoOfElements));
			printf("\n");
			af_print(moddims(FinalSphericalGrid(d_NoOfAnglesTheta / 2 - p, d_NoOfAnglesPhi / 2 + q, span), 1, d_NoOfElements));
			printf("\n");
			af_print(moddims(FinalSphericalGrid(d_NoOfAnglesTheta / 2 + p, d_NoOfAnglesPhi / 2 - q, span), 1, d_NoOfElements));
			printf("\n");
			af_print(moddims(FinalSphericalGrid(d_NoOfAnglesTheta / 2 + p, d_NoOfAnglesPhi / 2 + q, span), 1, d_NoOfElements));
			printf("\n");*/
			/**************************************************************************************************************************************************************************/

		}
		//printf("elapsed seconds: %g at iteration %d X-oriented Polar Slices\n", timer::stop(),q);

		if (q * 180 / d_NoOfAnglesPhi == 45)       // At 45 degrees its redundant to compute same in the next block  so skip it
			continue;
		//timer::start();

		for (int p = 1; p <= d_NoOfLevelsTheta; ++p)        // Y oriented pair of Polar slices, (d_NoOfAnglesPhi / 2 - q, d_NoOfAnglesPhi / 2 + q)
		{
			double angleTheta = p  * Pi / d_NoOfAnglesTheta;

			/**************************************************************************************************************************************************************************/
			// YY block
			alpha_factor = cos(angleTheta) * sin(anglePhi);				 // Scaling needed as defined in the paper
			beta_factor = cos(angleTheta) * cos(anglePhi);
			gamma_factor = sin(angleTheta);

			ReorderedImage_OperateColumns = Image3D;                                         // No swap directly we can operate on columns
			FrFT1D_Uniform_Image3D = ComputeNonVectorColumnwise_FrFTUniform(ReorderedImage_OperateColumns, beta_factor, d_NoOfElements); // This is the most expensive operation of order (N+1)^3 log(N+1)

			if (q == 1 && p == 1)         // Computing Polar Slice at anglephi = 0, only need to be computed once
			{ 
				array Central2D_YZSlice2D = FrFT1D_Uniform_Image3D(N / 2, span, span);
				Central2D_YZSlice2D = moddims(Central2D_YZSlice2D, d_NoOfElements, d_NoOfElements);
				//af_print(Central2D_YZSlice2D);

				array properOrientedSlice = flip(Central2D_YZSlice2D.T(), 0);          // Verified match  
				array Polar2D = Get2DFullPolarDFT(properOrientedSlice, d_NoOfElements, d_NoOfAnglesTheta, d_NoOfLevelsTheta);
				Polar2D = moddims(Polar2D, d_NoOfAnglesTheta, 1, d_NoOfElements);
				//af_print(Polar2D);
				FinalSphericalGrid(span, 0, span) = Polar2D;
				//af_print(BruteForceSphericalGrid(span, 0, span) - FinalSphericalGrid(span, 0, span));


				// Special Operation
				//af_print(Polar2D(d_NoOfAnglesTheta / 2, span));
				array SpecialLineZ = Polar2D(d_NoOfAnglesTheta / 2, span);
				//af_print(SpecialLineZ);
				array TiledZ = tile(SpecialLineZ, 1, d_NoOfAnglesPhi, 1);
				//af_print(TiledZ);
				//af_print(FinalSphericalGrid(d_NoOfAnglesTheta / 2, span, span));
				FinalSphericalGrid(d_NoOfAnglesTheta / 2, span, span) = TiledZ;        // This line is common to all grids !!! VERY VERY Special
			}

			ComputeNonVector_NonUniformFrFT_SingleColumn(FrFT1D_Uniform_Image3D.T(), alpha_factor, d_NoOfElements, Silce2D, Slice2D_Conj); //  2nd Swap, X-axis as columns now

			// Gather 4 lines
			ComputeFinal_2ComplementaryLines(Silce2D.T(), gamma_factor, d_NoOfElements, Line1, Line1_Conj);          // 3rd Swap , Z axis as columns
			ComputeFinal_2ComplementaryLines(Slice2D_Conj.T(), gamma_factor, d_NoOfElements, Line2, Line2_Conj);     // Swap Lines ! // 3rd Swap , Z axis as columns


			FinalSphericalGrid(p, d_NoOfAnglesPhi / 2 - q, span) = Line1;
			FinalSphericalGrid(p, d_NoOfAnglesPhi / 2 + q, span) = Line2;

			FinalSphericalGrid(d_NoOfAnglesTheta - p, d_NoOfAnglesPhi / 2 - q, span) = conjg(Line1_Conj);
			FinalSphericalGrid(d_NoOfAnglesTheta - p, d_NoOfAnglesPhi / 2 + q, span) = conjg(Line2_Conj);

			/*printf("\n My Special Solution 1st YY block\n");
			af_print(moddims(FinalSphericalGrid(p, d_NoOfAnglesPhi / 2 - q, span), 1, d_NoOfElements));
			printf("\n");
			af_print(moddims(FinalSphericalGrid(p, d_NoOfAnglesPhi / 2 + q, span), 1, d_NoOfElements));
			printf("\n");
			af_print(moddims(FinalSphericalGrid(d_NoOfAnglesTheta - p, d_NoOfAnglesPhi / 2 - q, span), 1, d_NoOfElements));
			printf("\n");
			af_print(moddims(FinalSphericalGrid(d_NoOfAnglesTheta - p, d_NoOfAnglesPhi / 2 + q, span), 1, d_NoOfElements));
			printf("\n");*/
			/**************************************************************************************************************************************************************************/

			//if (p * 180 / d_NoOfAnglesTheta == 45)       // At 45 degrees its redundant to compute in the next block       so skip it
			//	continue;

			
		}
		//printf("elapsed seconds: %g at iteration %d Y-oriented Polar Slices\n", timer::stop(),q);
	}
	
}

static int d_NoOfElements = 65;		/* IT is always odd */

std::vector<std::vector<std::string>> ReadMetaData()
{
	std::vector<std::vector<std::string>> result;
	std::ifstream myFile("../../../ModelNet40Preprocessing/archive/metadata_modelnet40.csv");

	if (!myFile.is_open()) throw std::runtime_error("Could not open file");

	std::string line, colstr;
	while (std::getline(myFile, line))
	{
		std::stringstream ss(line);
		std::vector<std::string> colVals;
		while (std::getline(ss, colstr, ',')) {
			colVals.push_back(colstr);
		}
		result.push_back(colVals);
	}
	myFile.close();
	return result;
}

void ReadSingleVolume(const std::string voxFileName, double* pReadVals)
{
	int count = 0;
	std::ifstream myFile(voxFileName);
	if (!myFile.is_open()) throw std::runtime_error("Could not open file");
	
	for (int x = 0; x < d_NoOfElements; x++) 
	{
		for (int y = 0; y < d_NoOfElements; y++)
		{
			for (int z = 0; z < d_NoOfElements; z++)
			{
				myFile >> pReadVals[count];
				count++;
			}
		}
	}
}

void FromVolumeComputeSphericalGrid(array& vol3D, array& voxSphericalGrid, size_t& numValues)
{
	// Initialization
	int B = (d_NoOfElements - 1) / 2;
	int d_NoOfAnglesTheta = 2 * B; /*Always even*/
	int d_NoOfAnglesPhi = B; /*Always even*/
	int d_NoOfLevelsTheta = ceil(float(d_NoOfAnglesTheta - 2) /
		4);  // This is where partitioning begins
	int d_NoOfLevelsPhi = ceil(float(d_NoOfAnglesPhi - 2) / 4);

	// Compute Spherical Grid
	voxSphericalGrid =
		constant(0, d_NoOfAnglesTheta, d_NoOfAnglesPhi, d_NoOfElements,
			PRECISION_COMPLEX);  // Angle Phi vs. Polar slices : Each
								 // Polar slice has No of elevation
								 // angle Theta vs. Radial data
	ComputeFastSphericalPolarFourierTransform(
		vol3D, voxSphericalGrid, d_NoOfElements, d_NoOfAnglesTheta,
		d_NoOfLevelsTheta, d_NoOfAnglesPhi, d_NoOfLevelsPhi);

	numValues = d_NoOfAnglesTheta * d_NoOfAnglesPhi * d_NoOfElements;
}

void WriteAbsSphericalGrid(const std::string absSpGridFileName, const size_t numVals, const double* d_Device)
{
	std::ofstream myFile(absSpGridFileName);
	for (int j = 0; j < numVals; ++j)
	{
		myFile << d_Device[j];
		myFile << "\n";
	}
	myFile.close();
}

// Using MATLAB these files were generated
void MatchExternalMATLABOutputAsGroundTruth()
{
	double* pReadVals = new double[d_NoOfElements * d_NoOfElements * d_NoOfElements];
	std::string voxFileName = "../../TestSerializedVolume.csv";

	// Load volume or voxel image
	ReadSingleVolume(voxFileName, pReadVals);
	array vol3D(d_NoOfElements, d_NoOfElements, d_NoOfElements, pReadVals);
	vol3D = vol3D.T();  // Needed as defined in the CUDAMexArrayFireSphericalPolarDFT file

	// Compute spherical polar fourier transform
	array voxSphericalGrid;
	size_t numValues;
	FromVolumeComputeSphericalGrid(vol3D, voxSphericalGrid, numValues);

	// ArrayFire computations to get final device pointer
	array voxSphericalGridReal = real(voxSphericalGrid);
	array voxSphericalGridImag = imag(voxSphericalGrid);
	array voxAbsFinalSphericalGrid = sqrt(pow(voxSphericalGridReal, 2) + pow(voxSphericalGridImag, 2));
	double* d_FinalAbsPolarGrid = voxAbsFinalSphericalGrid.host<double>();

	std::string spGridFileName = "../../TestSerializedSphericalGridArrayFire.csv";
	WriteAbsSphericalGrid(spGridFileName, numValues, d_FinalAbsPolarGrid);
	freeHost(d_FinalAbsPolarGrid);
}

using namespace concurrency;

int main(int argc, char* argv[]) 
{
	try {
		
		// This compares with MATLAB generated SPFT abs values from a random volume
		MatchExternalMATLABOutputAsGroundTruth();

		// Select a device and display arrayfire info
        int device = argc > 1 ? atoi(argv[1]) : 0;
        af::setDevice(device);
        af::info();

		// Load data
		std::vector<std::vector<std::string>> metaDataModelNet40 = ReadMetaData();
		
		std::vector<std::string> voxelsFolders;
		/*voxelsFolders.push_back("../../../ModelNet40Preprocessing/extracted_volumes/");
		voxelsFolders.push_back("../../../ModelNet40Preprocessing/extracted_volumes_noise/");*/
		/*voxelsFolders.push_back("../../../ModelNet40Preprocessing/extracted_volumes_missing/");
		voxelsFolders.push_back("../../../ModelNet40Preprocessing/extracted_volumes_outliers/");*/
		voxelsFolders.push_back("../../../ModelNet40Preprocessing/extracted_volumes_rot_noise/");

		double* pReadVals = new double[d_NoOfElements * d_NoOfElements * d_NoOfElements];

		// Output folder
		std::vector<std::string> sphericalGridFolders;
		/*sphericalGridFolders.push_back("../../../ModelNet40Preprocessing/extracted_spherical_grids/");
		sphericalGridFolders.push_back("../../../ModelNet40Preprocessing/extracted_spherical_grids_noise/");*/
		/*sphericalGridFolders.push_back("../../../ModelNet40Preprocessing/extracted_spherical_grids_missing/");
		sphericalGridFolders.push_back("../../../ModelNet40Preprocessing/extracted_spherical_grids_outliers/");*/
		sphericalGridFolders.push_back("../../../ModelNet40Preprocessing/extracted_spherical_grids_rot_noise/");
		
		for(int fIndex = 0; fIndex < sphericalGridFolders.size(); fIndex++)
			fs::create_directory(sphericalGridFolders[fIndex]);

		int rowOffset = 1;  // Remove this line once done
		for (int fIndex = 0; fIndex < sphericalGridFolders.size(); fIndex++)
		{
			//parallel_for(size_t(0), metaDataModelNet40.size(), [&](size_t row) - parallel doesn't work due to ArrayFire
			for (int row = 1; row < metaDataModelNet40.size(); row++) // Remove this and replace it with 
			{
				timer::start();
				std::string strObj = metaDataModelNet40[row][0];
				std::string strClass = metaDataModelNet40[row][1];
				std::string strSplit = metaDataModelNet40[row][2];

				std::string voxFileName = voxelsFolders[fIndex]
					+ strClass + "/" + strSplit + "/" + strObj + ".csv";

				// Load volume or voxel image
				ReadSingleVolume(voxFileName, pReadVals);
				array vol3D(d_NoOfElements, d_NoOfElements, d_NoOfElements, pReadVals);
				vol3D = vol3D.T();  // Needed as defined in the CUDAMexArrayFireSphericalPolarDFT file

				// Compute spherical polar fourier transform
				array voxSphericalGrid;
				size_t numValues;
				FromVolumeComputeSphericalGrid(vol3D, voxSphericalGrid, numValues);

				// ArrayFire computations to get final device pointer
				array voxSphericalGridReal = real(voxSphericalGrid);
				array voxSphericalGridImag = imag(voxSphericalGrid);
				array voxAbsFinalSphericalGrid = sqrt(pow(voxSphericalGridReal, 2) + pow(voxSphericalGridImag, 2));
				double* d_FinalAbsPolarGrid = voxAbsFinalSphericalGrid.host<double>();

				// Save in right directory
				std::string spGridClassFolder = sphericalGridFolders[fIndex] + strClass + "/";
				fs::create_directory(spGridClassFolder);
				std::string spGridInnerSplitFolder = sphericalGridFolders[fIndex] + strClass + "/" + strSplit + "/";
				fs::create_directory(spGridInnerSplitFolder);

				// Write csv
				std::string spGridFileName = sphericalGridFolders[fIndex]
					+ strClass + "/" + strSplit + "/" + strObj + "_sp.csv";
				WriteAbsSphericalGrid(spGridFileName, numValues, d_FinalAbsPolarGrid);
				freeHost(d_FinalAbsPolarGrid);
				printf("%d, elapsed seconds: %g\n", row, timer::stop());
			}
			//});
		}

		delete[] pReadVals;
	}
	catch (af::exception& e) {
		fprintf(stderr, "%s\n", e.what());
		throw;
	}
	return 0;
}