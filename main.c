/*
 * main.c
 *
 *  Created on: Jul 10, 2012
 *      Author: toto
 */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <fftw3.h>

//fft and ifft for float complex number
fftwf_complex *fftwf_data(fftwf_complex *input, int ndata, int nfft);
fftwf_complex *ifftwf_data(fftwf_complex *fdata, int ndata, int nfft);

//fft and ifft for double complex number
fftw_complex *fftwd_data(fftw_complex *input, int ndata, int nfft);
fftw_complex *ifftwd_data(fftw_complex *fdata, int ndata, int nfft);
int main( int argc, char** argv )
{
	fftwf_complex *data;
	fftwf_complex *fdata;
	fftwf_complex *idata;
	int ndata, nfft, i;
	unsigned int seed = 1234;

	if(argc<2)
		return(1);
	sscanf(argv[1], "%i", &ndata);
	sscanf(argv[2], "%i", &nfft);

	//allocate memory for data
	data = ( fftwf_complex* ) fftwf_malloc( sizeof( fftwf_complex ) * ndata );

	// CREATE INPUT DATA
	srand ( seed );
	for( i=0; i<ndata; i++ ) {
		data[i][0] = (double) (rand()/12345678);	//real data
		data[i][1] = 0.0;		//imaginer data
	}
	//print input data
	printf("Input Data \n");
	for( i = 0 ; i < ndata ; i++ )
		printf( "%i.  %2.5f  %2.5f \n", i, data[i][0], data[i][1]);

	// FFT PROCESS
	fdata = fftwf_data(data, ndata, nfft);
	printf("\nFFT Result \n");
	for( i = 0 ; i < nfft ; i++ )
		printf( "%i.  %2.5f  %2.5f \n", i, fdata[i][0], fdata[i][1]);

	// INVERS FFT PROCESS
	idata = ifftwf_data(fdata, ndata, nfft);
	printf("\nIFFT Result \n");
	for( i = 0 ; i < ndata ; i++ )
		printf( "%i.  %2.5f  %2.5f \n", i, idata[i][0], idata[i][1]);
	printf("\n");

	fftwf_free( data );
	fftwf_free( fdata );
	fftwf_free( idata );

	return 0;
}

/* --------------------------------------------------------------------
 * FFT for floating complex number data
 * Comparation with Matlab Command :
 * b = fft(A) --> b = fftwf_data(A, lenA, lenA)
 * b = fft(A, nfft) --> b = fftwf_data(A, lenA, nfft)
 * -------------------------------------------------------------------- */
fftwf_complex *fftwf_data(fftwf_complex *input, int ndata, int nfft)
{
	fftwf_complex *paddata;
	fftwf_complex *fdata;
	fftwf_plan plan_forward;

	if(nfft<ndata) {
		fprintf(stderr, "nfft < ndata \n");
		exit(0);
	}

	//allocate memory for data+padding
	paddata = ( fftwf_complex* ) fftwf_malloc( sizeof( fftwf_complex ) * nfft );

	//allocate data for output fft process
	fdata = ( fftwf_complex* ) fftwf_malloc( sizeof( fftwf_complex ) * nfft );

	//padding input data
	memset(paddata[0], 0, nfft*sizeof(fftwf_complex));
	memcpy(paddata[0], input[0], ndata*sizeof(fftwf_complex));

	plan_forward  = fftwf_plan_dft_1d( nfft, paddata, fdata, FFTW_FORWARD, FFTW_ESTIMATE );
	fftwf_execute( plan_forward ); //execute fft

	//free allocate memory
	fftwf_destroy_plan( plan_forward );
	fftwf_free( paddata );

	return(fdata);
}

/* --------------------------------------------------------------------
 * IFFT for floating complex number data
 * Comparation with Matlab Command :
 * b = ifft(A) --> b = ifftwf_data(A, lenA, lenA)
 * b = ifft(A, nfft) --> b = ifftwf_data(A, lenA, nfft)
 * -------------------------------------------------------------------- */
fftwf_complex *ifftwf_data(fftwf_complex *fdata, int ndata, int nfft)
{
	fftwf_complex *idata;
	fftwf_complex *idatapad;
	fftwf_plan plan_backward;
	int i;

	if(nfft<ndata) {
		fprintf(stderr, "nfft < ndata \n");
		exit(0);
	}

	//allocate memory for data
	idata = ( fftwf_complex* ) fftwf_malloc( sizeof( fftwf_complex ) * ndata );

	//allocate memory for inversfft + padding
	idatapad = ( fftwf_complex* ) fftwf_malloc( sizeof( fftwf_complex ) * nfft );

	plan_backward = fftwf_plan_dft_1d( nfft, fdata, idatapad, FFTW_BACKWARD, FFTW_ESTIMATE );
	fftwf_execute( plan_backward ); //execute invers fft

	//get invers fft with length ndata and divide the value with nfft
	for( i=0; i<ndata; i++ ){
		idata[i][0] = idatapad[i][0]/nfft;	//real data
		idata[i][1] = idatapad[i][1]/nfft;	//imaginer data
	}

	//free allocate memory
	fftwf_destroy_plan( plan_backward );
	fftwf_free( idatapad );

	return(idata);
}

/* --------------------------------------------------------------------
 * FFT for double complex number data
 * Comparation with Matlab Command :
 * b = fft(A) --> b = fftwf_data(A, lenA, lenA)
 * b = fft(A, nfft) --> b = fftwd_data(A, lenA, nfft)
 * -------------------------------------------------------------------- */
fftw_complex *fftwd_data(fftw_complex *input, int ndata, int nfft)
{
	fftw_complex *paddata;
	fftw_complex *fdata;
	fftw_plan plan_forward;

	if(nfft<ndata) {
		fprintf(stderr, "nfft < ndata \n");
		exit(0);
	}

	//allocate memory for data+padding
	paddata = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * nfft );

	//allocate data for output fft process
	fdata = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * nfft );

	//padding input data
	memset(paddata[0], 0, nfft*sizeof(fftw_complex));
	memcpy(paddata[0], input[0], ndata*sizeof(fftw_complex));

	plan_forward  = fftw_plan_dft_1d( nfft, paddata, fdata, FFTW_FORWARD, FFTW_ESTIMATE );
	fftw_execute( plan_forward ); //execute fft

	//free allocate memory
	fftw_destroy_plan( plan_forward );
	fftw_free( paddata );

	return(fdata);
}

/* --------------------------------------------------------------------
 * IFFT for double complex number data
 * Comparation with Matlab Command :
 * b = ifft(A) --> b = ifftwd_data(A, lenA, lenA)
 * b = ifft(A, nfft) --> b = ifftwd_data(A, lenA, nfft)
 * -------------------------------------------------------------------- */
fftw_complex *ifftwd_data(fftw_complex *fdata, int ndata, int nfft)
{
	fftw_complex *idata;
	fftw_complex *idatapad;
	fftw_plan plan_backward;
	int i;

	if(nfft<ndata) {
		fprintf(stderr, "nfft < ndata \n");
		exit(0);
	}

	//allocate memory for data
	idata = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * ndata );

	//allocate memory for inversfft + padding
	idatapad = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * nfft );

	plan_backward = fftw_plan_dft_1d( nfft, fdata, idatapad, FFTW_BACKWARD, FFTW_ESTIMATE );
	fftw_execute( plan_backward ); //execute invers fft

	//get invers fft with length ndata and divide the value with nfft
	for( i=0; i<ndata; i++ ){
		idata[i][0] = idatapad[i][0]/nfft;	//real data
		idata[i][1] = idatapad[i][1]/nfft;	//imaginer data
	}

	//free allocate memory
	fftw_destroy_plan( plan_backward );
	fftw_free( idatapad );

	return(idata);
}

