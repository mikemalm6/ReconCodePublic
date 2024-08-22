/*	Gridding Routines.  To be included as header files. */
/* Michael Malmberg edited from Brian Hargreaves gridkb to split up the 
precomputation of weights from the actual gridding itself */

#include <math.h>
#include <stdio.h>

void gridprep(kx,ky,nsamples,
	indxminall, indxmaxall, indyminall, indymaxall, wts, gridsize, convwidth, kerneltable, nkernelpts,nwt)

/* Gridding function that uses a lookup table for a circularly
	symmetric convolution kernel, with linear interpolation.  
	See Notes below */

/* INPUT/OUTPUT */

double *kx;		/* Array of kx locations of samples. */
double *ky;		/* Array of ky locations of samples. */
int nsamples;		/* Number of k-space samples, total. */
double *indxminall;
double *indxmaxall;
double *indyminall;
double *indymaxall;
double *wts;
int nwt;
int gridsize;		/* Number of points in kx and ky in grid. */
double convwidth;	/* Kernel width, in grid points.	*/
double *kerneltable;	/* 1D array of convolution kernel values, starting
				at 0, and going to convwidth. */
int nkernelpts;		/* Number of points in kernel lookup-table */

{
int kcount;		/* Counts through k-space sample locations */
int gcount1, gcount2;	/* Counters for loops */
long int wcount;
int col;		/* Grid Columns, for faster lookup. */

double kwidth;			/* Conv kernel width, in k-space units */
double dkx,dky,dk;		/* Delta in x, y and abs(k) for kernel calc.*/
int ixmin,ixmax,iymin,iymax;	/* min/max indices that current k may affect*/
int kernelind;			/* Index of kernel value, for lookup. */
double fracdk;			/* Fractional part of lookup index. */
double fracind;			/* Fractional lookup index. */
double kern;			/* Kernel value, avoid duplicate calculation.*/
double *sgrptr;			/* Aux. pointer, for loop. */
double *sgiptr;			/* Aux. pointer, for loop. */
double *indxminaptr;
double *indxmaxaptr;
double *indyminaptr;
double *indymaxaptr;
double *wtsptr;
double nwtd = (double)nwt;
double nsampleswts = nwt*nsamples;
int gridsizesq;			/* Square of gridsize */
int gridcenter;			/* Index in output of kx,ky=0 points. */
int gptr_cinc;			/* Increment for grid pointer. */
int wtsctr;
wtsctr = 0;
gridcenter = gridsize/2;	/* Index of center of grid. */
kwidth = convwidth/(double)(gridsize);	/* Width of kernel, in k-space units. */


/* ========= Zero Output Points ========== */

indxminaptr = indxminall;
indxmaxaptr = indxmaxall;
indyminaptr = indyminall;
indymaxaptr = indymaxall;
wtsptr = wts;
gridsizesq = gridsize*gridsize;
 
for (kcount = 0; kcount < nsamples; kcount++)
{
    *indxminaptr++ = 0.0;
    *indxmaxaptr++ = 0.0;
    *indyminaptr++ = 0.0;
    *indymaxaptr++ = 0.0;
}
for (wcount = 0; wcount < nsampleswts; wcount++) {
    *wtsptr++ = 0.0;
}
indxminaptr = indxminall;
indxmaxaptr = indxmaxall;
indyminaptr = indyminall;
indymaxaptr = indymaxall;
wtsptr = wts;
/* ========= Loop Through k-space Samples ========= */
                
for (kcount = 0; kcount < nsamples; kcount++)
	{

	/* ----- Find limit indices of grid points that current
		 sample may affect (and check they are within grid) ----- */

	ixmin = (int) ((*kx-kwidth)*gridsize +gridcenter);
	if (ixmin < 0) ixmin=0;
	ixmax = (int) ((*kx+kwidth)*gridsize +gridcenter)+1;
	if (ixmax >= gridsize) ixmax=gridsize-1;
	iymin = (int) ((*ky-kwidth)*gridsize +gridcenter);
	if (iymin < 0) iymin=0;
	iymax = (int) ((*ky+kwidth)*gridsize +gridcenter)+1;
	if (iymax >= gridsize) iymax=gridsize-1;
    
    *indxminaptr = (double)ixmin;
    *indxmaxaptr = (double)ixmax;
    *indyminaptr = (double)iymin;
    *indymaxaptr = (double)iymax;

		
	  /* Increment for grid pointer at end of column to top of next col.*/
	gptr_cinc = gridsize-(iymax-iymin)-1;	/* 1 b/c at least 1 sgrptr++ */
		
	wtsctr = 0;
	for (gcount1 = ixmin; gcount1 <= ixmax; gcount1++)
        	{
		dkx = (double)(gcount1-gridcenter) / (double)gridsize  - *kx;
		for (gcount2 = iymin; gcount2 <= iymax; gcount2++)
			{
            dky = (double)(gcount2-gridcenter) / 
            (double)gridsize - *ky;

			dk = sqrt(dkx*dkx+dky*dky);	/* k-space separation*/

			if (dk < kwidth)	/* sample affects this grid point	*/
			    {
				/* Find index in kernel lookup table */
			    fracind = dk/kwidth*(double)(nkernelpts-1);
			    kernelind = (int)fracind;
			    fracdk = fracind-(double)kernelind;

				/* Linearly interpolate in kernel lut */
			    kern = kerneltable[(int)kernelind]*(1-fracdk)+
			    		kerneltable[(int)kernelind+1]*fracdk;

                *wtsptr = kern;
                
			    }
            wtsctr += 1;
            wtsptr++;
			}
		}
	kx++;		/* Advance kx pointer */
	ky++;   	/* Advance ky pointer */
    indxminaptr++;
    indxmaxaptr++;
    indyminaptr++;
    indymaxaptr++;
    wtsptr = wtsptr + (nwt-wtsctr);
	}
}


