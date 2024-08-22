
#include "mex.h" 
#include "gridfastroutine.c"
#include "matrix.h"
#include "math.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])

{
double *s_real;		/* real-part of data samples.	*/
double *s_imag;		/* imaginary-part of data samples. */
int nsamples;		/* number of data samples (or trajectory points) */
double *dens;		/* (pre) density correction factor at each traj. loc.*/
double *sg_real;	/* real-part of gridded data sample. */
double *sg_imag;	/* imaginary-part of gridded data sample. */
int gridsize;		/* size of grid, in samples.  ie gridsize x gridsize */
double *ixminall;
double *ixmaxall;
double *iyminall;
double *iymaxall;
double *wtsall;

mxArray *dens_mat;	
nsamples = mxGetM(prhs[0]) * mxGetN(prhs[0]);	/* Samples may be passed as
							1xN, Nx1 or 2D array */

s_real = mxGetPr(prhs[0]);	/* Get real parts of data samples. */
s_imag = mxGetPi(prhs[0]);	/* Get imaginary parts of data samples. */
dens = mxGetPr(prhs[1]);	/* Get density correction factors. */
gridsize = (int)(*mxGetPr(prhs[2]));	/* Get grid size. */
ixminall = mxGetPr(prhs[3]);
ixmaxall = mxGetPr(prhs[4]);
iyminall = mxGetPr(prhs[5]);
iymaxall = mxGetPr(prhs[6]);
wtsall = mxGetPr(prhs[7]);

int nwt = (int)mxGetM(prhs[7]);

plhs[0] = mxCreateDoubleMatrix(gridsize,gridsize,mxCOMPLEX);	/* output */              

sg_real = mxGetPr(plhs[0]);
sg_imag = mxGetPi(plhs[0]);

gridfast(s_real,s_imag,nsamples, dens,
		sg_real,sg_imag, ixminall, ixmaxall, iyminall, iymaxall,wtsall,gridsize,nwt);

}





