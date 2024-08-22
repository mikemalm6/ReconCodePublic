
#include "mex.h" 
#include "gridpreproutine.c"
#include "matrix.h"
#include "math.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])


{
double *kx;		/* k-space trajectory x locations.	*/
double *ky;		/* k-space trajectory y locations.	*/

int nsamples;		/* number of data samples (or trajectory points) */

int gridsize;		/* size of grid, in samples.  ie gridsize x gridsize */
double convwidth;	/* width of convolution kernel for gridding. */
double *kerneltable;	/* lookup-table for linearly-interpolated kernel.  */
int nkernelpts;		/* number of points in kernel lookup table. */
double *indxminall;
double *indxmaxall;
double *indyminall;
double *indymaxall;
double *wts;
mwSize nR;
mwSize nC;
mxArray *dens_mat;
nR = mxGetM(prhs[0]);
nC = mxGetN(prhs[0]);
nsamples = nR * nC;	/* Samples may be passed as
							1xN, Nx1 or 2D array */
kx = mxGetPr(prhs[0]);		/* Get kx locations */
ky = mxGetPi(prhs[0]);		/* Get ky locations */
gridsize = (int)(*mxGetPr(prhs[1]));	/* Get grid size. */
	
if (nrhs > 2)
	convwidth = *mxGetPr(prhs[2]);			/* Get conv width */
else
	convwidth = 1.0;				/* Assign default. */

int nwt = (int)ceil((convwidth+1)*(convwidth+1)*4*1.125);
kerneltable = mxGetPr(prhs[3]);				/* Get kernel table.*/
nkernelpts = mxGetM(prhs[3]) * mxGetN(prhs[3]);		/* and # points. */

plhs[0] = mxCreateDoubleMatrix(nR,nC,mxREAL);
if (nlhs >= 2) {
    plhs[1] = mxCreateDoubleMatrix(nR,nC,mxREAL);
    if (nlhs >=3) {
        plhs[2] = mxCreateDoubleMatrix(nR,nC,mxREAL);
        if (nlhs >=4) {
            plhs[3] = mxCreateDoubleMatrix(nR,nC,mxREAL);
            if (nlhs >=5) {
                plhs[4] = mxCreateDoubleMatrix(nwt,nsamples,mxREAL); 
            }
        }
    }
}
                
indxminall = mxGetPr(plhs[0]);
indxmaxall = mxGetPr(plhs[1]);
indyminall = mxGetPr(plhs[2]);
indymaxall = mxGetPr(plhs[3]);
wts = mxGetPr(plhs[4]);

gridprep(kx,ky,nsamples,
		indxminall, indxmaxall, indyminall, indymaxall,wts,gridsize, convwidth, kerneltable, nkernelpts,nwt);

}





