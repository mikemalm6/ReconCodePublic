#include <math.h>
#include <stdio.h>

void gridfast(s_real,s_imag,nsamples, dcf,
	sg_real,sg_imag, ixminall, ixmaxall, iyminall, iymaxall, wtsall, gridsize,nwt)

/* Gridding function that uses a lookup table for a circularly
	symmetric convolution kernel, with linear interpolation.  
	See Notes below */

/* INPUT/OUTPUT */

double *s_real;		/* Sampled data, real part. */
double *s_imag; 	/* Sampled data, real part. */
int nsamples;		/* Number of k-space samples, total. */
double *dcf;		/* Density compensation factors. */

double *sg_real;	/* OUTPUT array, real parts of data. */
double *sg_imag;	/* OUTPUT array, imag parts of data. */	 
double *ixminall;
double *ixmaxall;
double *iyminall;
double *iymaxall;
double *wtsall;

int gridsize;		/* Number of points in kx and ky in grid. */
int nwt;


/*------------------------------------------------------------------
	NOTES:

	This uses the following formula, which describes the contribution
	of each data sample to the value at each grid point:

		grid-point value += data value * dcf * kernel(dk)

	where:
		data value is the complex sampled data value.
		dcf is the density compensation factor for the sample point.
		kernel is the convolution kernel function.
		dk is the k-space separation between sample point and
			grid point.

	"grid units"  are integers from 0 to gridsize-1, where
	the grid represents a k-space of -.5 to .5.

	The convolution kernel is passed as a series of values 
	that correspond to "kernel indices" 0 to nkernelpoints-1.
  ------------------------------------------------------------------ */

{
int kcount;		/* Counts through k-space sample locations */
int gcount1, gcount2;	/* Counters for loops */




int ixmin,ixmax,iymin,iymax;	/* min/max indices that current k may affect*/
double kern;			/* Kernel value, avoid duplicate calculation.*/
double *sgrptr;			/* Aux. pointer, for loop. */
double *sgiptr;			/* Aux. pointer, for loop. */

int gridsizesq;			/* Square of gridsize */
int gridcenter;			/* Index in output of kx,ky=0 points. */
int gptr_cinc;			/* Increment for grid pointer. */
int wtsctr;
int offset;
wtsctr = 0;
gridcenter = gridsize/2;	/* Index of center of grid. */
double srdcf = 0.0;
double sidcf = 0.0;


/* ========= Zero Output Points ========== */

sgrptr = sg_real;
sgiptr = sg_imag;
gridsizesq = gridsize*gridsize;
 
for (gcount1 = 0; gcount1 < gridsizesq; ++gcount1)
        {
	*sgrptr++ = 0;
	*sgiptr++ = 0;
	}

/* ========= Loop Through k-space Samples ========= */
                
for (kcount = 0; kcount < nsamples; ++kcount)
	{

	/* ----- Find limit indices of grid points that current
		 sample may affect (and check they are within grid) ----- */

    
    ixmin = (int)*ixminall;
    ixmax = (int)*ixmaxall;
    iymin = (int)*iyminall;
    iymax = (int)*iymaxall;
	
	  /* Increment for grid pointer at end of column to top of next col.*/
	gptr_cinc = gridsize-(iymax-iymin)-1;	/* 1 b/c at least 1 sgrptr++ */
    offset = ixmin*gridsize+iymin;
	sgrptr = sg_real + offset;
	sgiptr = sg_imag + offset;
    srdcf = *s_real * *dcf;
    sidcf = *s_imag * *dcf;
	wtsctr = 0;
	for (gcount1 = ixmin; gcount1 <= ixmax; ++gcount1)
        {
		for (gcount2 = iymin; gcount2 <= iymax; ++gcount2)
            {
            kern = *wtsall;
            if (kern > 0) {
                *sgrptr += kern * srdcf;
                *sgiptr += kern * sidcf;
            }
            wtsctr += 1;
            ++wtsall;
			++sgrptr;
			++sgiptr;
        }
		sgrptr+= gptr_cinc;
		sgiptr+= gptr_cinc;
    }

	++dcf;		/* Advance dcf pointer */		
	++s_real;	/* Advance real-sample pointer */
	++s_imag;	/* Advance imag-sample pointer */    
    ++ixminall;
    ++ixmaxall;
    ++iyminall;
    ++iymaxall;
    wtsall = wtsall + (nwt-wtsctr);
	}
}
